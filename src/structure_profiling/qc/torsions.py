'''
Chain torsion / backbone-angle analysis for one Bio.PDB chain.

Collected here from src/structure_profiling (previously spread over
geometry/geometry_angle.py, used by geometry/geometry.py and
qc/qc_screen.py):
- dihedral_angle / bond_angle:  point-based torsion / angle utilities
- evaluate_ramachandran:        single-box phi/psi classification
- bond_angles:                  raw backbone bond angles (N-CA-C, CA-C-N;
                                ideal ~111 / ~117 degrees)
- torsion_angles:               phi / psi / omega per residue with
                                cis/trans conformation and Ramachandran
                                label (cis when |omega| < 90 deg)
- bond_angle_z:                 Z-scores of backbone bond angles vs the
                                Engh & Huber ideals (|Z| > 5 outlier)
- peptide_planarity:            wwPDB planarity per peptide unit
                                (omega dev > 30 deg, out-of-plane > 10 deg)
- phi_psi_series / rama_class / residue-type aware Ramachandran
  rama_molprobity:              (MolProbity style: Gly / Pro / pre-Pro)
- rama_grid_summary:            empirical-grid Ramachandran (RamaGrid)

run() returns the summary exported by pipelines/qc/chain_torsions.py to
the SQL table chain_torsions; check() gates it on the Torsions cutoffs
(the BACKBONE_GATES / STRUCTURE_GATES values of the former qc_screen
screens and the torsion gate in Geometry.check()).

Cutoffs:
    Z_OUTLIER_CUT = 5.0        |Z| above this flags a bond-angle outlier
    CIS_CUT = 90 deg           |omega| below this is a cis peptide
    OMEGA_DEV_CUT = 30 deg     deviation of omega from the nearest
                               cis(0) / trans(180) flags a peptide
                               planarity outlier
    MAINCHAIN_OOP_CUT = 10 deg out-of-plane N(i+1) flags a main-chain
                               planarity outlier
    GATES: angle_rmsz = 1.5, z_outlier_frac = 0.01,
    torsion_null_frac = 0.05, peptide_outlier_frac = 0.05,
    mainchain_outlier_frac = 0.05, rama_outlier_frac = 0.05

references:
wwPDB OneDep validation reports (X-ray + EM guides, updated 2025);
Engh & Huber bond-angle targets; MolProbity (Chen et al. 2010,
Read et al. 2011).
'''
import numpy as np
import pandas as pd
from Bio import PDB

from ..geometry.geometry_angle import GeometryAngle


class Torsions:
    '''backbone torsion / bond-angle analysis of one Bio.PDB chain.'''

    # ideal backbone bond angles (mean, sd) in degrees, Engh & Huber style -
    # canonical values live in GeometryAngle (geometry/geometry_angle.py)
    IDEAL_ANGLES = GeometryAngle.IDEAL_ANGLES
    Z_OUTLIER_CUT: float = GeometryAngle.Z_OUTLIER_CUT
    CIS_CUT: float = 90.0            # deg, |omega| below this is cis
    OMEGA_DEV_CUT: float = 30.0      # deg, omega dev above this is outlier
    MAINCHAIN_OOP_CUT: float = 10.0  # deg, out-of-plane above this fails

    # gates (fail -> reason), the former BACKBONE_GATES /
    # STRUCTURE_GATES of qc_screen.py and the torsion gate in Geometry.check()
    GATES = {
        'angle_rmsz': 1.5,              # RMS |Z| over backbone bond angles
        'z_outlier_frac': 0.01,         # fraction of |Z|>5 angle outliers
        'torsion_null_frac': 0.05,      # residues without measurable phi/psi
        'peptide_outlier_frac': 0.05,   # omega dev > OMEGA_DEV_CUT
        'mainchain_outlier_frac': 0.05, # out-of-plane > MAINCHAIN_OOP_CUT
        'rama_outlier_frac': 0.05,      # MolProbity-style Rama outliers
    }

    def __init__(self, chain):
        self.chain = chain

    # ------------------------------------------------------------------
    # residue selection (same filter as Breaks / GeometryDistance)
    # ------------------------------------------------------------------

    def standard_residues(self):
        '''standard amino-acid residues (skip hetero / water).

        Cutoff: none - filters on residue id[0] == ' ' only.
        '''
        return [r for r in self.chain if r.id[0] == " "]

    # ------------------------------------------------------------------
    # point-based angle / torsion utilities
    # ------------------------------------------------------------------

    @staticmethod
    def dihedral_angle(p0, p1, p2, p3):
        '''dihedral angle (degrees) between four points p0-p1-p2-p3.

        Cutoff: none - raw signed angle in (-180, 180].
        '''
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

    @staticmethod
    def bond_angle(a, b, c):
        '''bond angle (degrees) a-b-c around the central atom b.

        Cutoff: none - raw angle in [0, 180].
        '''
        v1, v2 = a - b, c - b
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))

    @staticmethod
    def evaluate_ramachandran(phi, psi):
        '''
        Check if phi/psi angles are in allowed regions
        General allowed regions for L-amino acids (degrees)
        https://jenalib.leibniz-fli.de/ImgLibDoc/glossary/IMAGE_RAMA.htm

        Cutoffs (analytic boxes, degrees): beta-sheet phi in (-180, -60)
        and psi in (90, 180); right-handed helix phi in (-120, -30), psi
        in (-60, 30); left-handed helix phi in (45, 90), psi in (0, 60);
        'allowed' phi in [-180, 0], psi in [0, 180]; else 'bad'.
        phi/psi are expected in degrees; None / NaN -> 'null'
        '''
        if phi is None or psi is None:
            return 'null'
        if not (np.isfinite(phi) and np.isfinite(psi)):
            return 'null'

        # Core allowed regions (approximate)
        # Beta-sheet: parallel or antiparallel
        if -180 < phi < -60 and 90 < psi < 180:
            return "beta-sheet"
        # Alpha-helix region
        elif -120 < phi < -30 and -60 < psi < 30:
            return "right-handed-helix"
        # Left-handed helix (rare, usually only Gly)
        elif 45 < phi < 90 and 0 < psi < 60:
            return "left-handed-helix"
        # Allowed but not favorable
        elif -180 <= phi <= 0 and 0 <= psi <= 180:
            return "allowed"
        return "bad"

    # ------------------------------------------------------------------
    # backbone bond angles
    # ------------------------------------------------------------------

    def bond_angles(self):
        '''
        backbone bond angles (degrees): N-CA-C within each residue and
        CA-C-N across the peptide bond; ideal values ~111 / ~117 degrees.
        Cutoff: none here - raw angles (bond_angle_z() flags |Z| > 5.0).
        returns DataFrame(residue, res_pos, angle_N_CA_C, angle_CA_C_N);
        angle_CA_C_N is None at peptide fragment ends.
        '''
        ppb = PDB.CaPPBuilder()
        res = []
        for pp in ppb.build_peptides(self.chain, aa_only=False):
            for i, residue in enumerate(pp):
                row = {
                    'residue': residue.resname,
                    'res_pos': residue.id[1],
                    'angle_N_CA_C': None,
                    'angle_CA_C_N': None,
                }
                if all(a in residue for a in ("N", "CA", "C")):
                    row['angle_N_CA_C'] = self.bond_angle(
                        residue["N"].coord, residue["CA"].coord, residue["C"].coord
                    )
                if i < len(pp) - 1 and all(a in residue for a in ("CA", "C")) \
                        and "N" in pp[i+1]:
                    row['angle_CA_C_N'] = self.bond_angle(
                        residue["CA"].coord, residue["C"].coord, pp[i+1]["N"].coord
                    )
                res.append(row)
        return pd.DataFrame(res)

    def bond_angle_z(self) -> pd.DataFrame:
        '''Z-scores of backbone bond angles vs IDEAL_ANGLES (|Z| > 5 outlier).

        Cutoff: Z_OUTLIER_CUT=5.0 on |Z|; angle RMS-Z expected in [0, 1]
        (wwPDB), former BACKBONE_GATES tolerates angle_rmsz up to 1.5.
        Covers intra-residue N-CA-C plus CA(i)-C(i)-N(i+1) and
        C(i)-N(i+1)-CA(i+1) bridge angles. Returns columns res_pos, resname,
        atoms, observed, ideal, sd, z, outlier.
        '''
        cols = ['res_pos', 'resname', 'atoms', 'observed',
                'ideal', 'sd', 'z', 'outlier']
        rows = []
        res = self.standard_residues()
        for i, r in enumerate(res):
            if all(a in r for a in ('N', 'CA', 'C')):
                ang = self.bond_angle(r['N'].coord, r['CA'].coord,
                                      r['C'].coord)
                mean, sd = self.IDEAL_ANGLES[('N', 'CA', 'C')]
                z = (ang - mean) / sd
                rows.append([r.id[1], r.resname, 'N-CA-C', ang, mean, sd,
                             z, abs(z) > self.Z_OUTLIER_CUT])
            if i < len(res) - 1:
                nxt = res[i + 1]
                if all(a in r for a in ('CA', 'C')) \
                        and all(a in nxt for a in ('N', 'CA')):
                    a1 = self.bond_angle(r['CA'].coord, r['C'].coord,
                                         nxt['N'].coord)
                    m1, s1 = self.IDEAL_ANGLES[('CA', 'C', 'N')]
                    rows.append([r.id[1], r.resname, 'CA-C-N', a1, m1, s1,
                                 (a1 - m1) / s1,
                                 abs((a1 - m1) / s1) > self.Z_OUTLIER_CUT])
                    a2 = self.bond_angle(r['C'].coord, nxt['N'].coord,
                                         nxt['CA'].coord)
                    m2, s2 = self.IDEAL_ANGLES[('C', 'N', 'CA')]
                    rows.append([r.id[1], r.resname, 'C-N-CA', a2, m2, s2,
                                 (a2 - m2) / s2,
                                 abs((a2 - m2) / s2) > self.Z_OUTLIER_CUT])
        return pd.DataFrame(rows, columns=cols)

    # ------------------------------------------------------------------
    # backbone torsion angles (phi / psi / omega)
    # ------------------------------------------------------------------

    def torsion_angles(self):
        '''backbone torsion angles (degrees)
        Calculate: φ (phi), ψ (psi), ω (omega) and look for Ramachandran
        outliers. omega is None at the end of each peptide fragment;
        conformation is cis (|omega| < CIS_CUT deg) or trans.
        Cutoff: cis when |omega| < 90 deg, else trans (pro cis is common,
        non-Pro cis is rare and counted in run() omega stats).
        '''
        ppb = PDB.CaPPBuilder()
        res = []
        for pp in ppb.build_peptides(self.chain, aa_only=False):
            phi_psi = pp.get_phi_psi_list()
            for i, (residue, (phi, psi)) in enumerate(zip(pp, phi_psi)):
                omega, conformation = None, 'null'
                if i < len(pp) - 1:
                    nxt = pp[i + 1]
                    if all(a in residue for a in ("CA", "C")) \
                            and all(a in nxt for a in ("N", "CA")):
                        omega = self.dihedral_angle(
                            residue["CA"].coord, residue["C"].coord,
                            nxt["N"].coord, nxt["CA"].coord,
                        )
                        conformation = 'cis' if abs(omega) < self.CIS_CUT \
                            else 'trans'
                # get_phi_psi_list returns radians; classify in degrees
                phi_d = None if phi is None else np.degrees(phi)
                psi_d = None if psi is None else np.degrees(psi)
                ramachandran = self.evaluate_ramachandran(phi_d, psi_d)
                res.append({
                    'residue': residue.resname,
                    'res_pos': residue.id[1],
                    'phi': phi_d,
                    'psi': psi_d,
                    'omega': omega,
                    'conformation': conformation,
                    'ramachandran': ramachandran,
                })
        return pd.DataFrame(res)

    # ------------------------------------------------------------------
    # peptide / main-chain planarity (wwPDB)
    # ------------------------------------------------------------------

    def peptide_planarity(self) -> pd.DataFrame:
        '''wwPDB planarity per peptide unit (i -> i+1).

        Cutoffs (wwPDB): peptide_outlier when omega_dev > OMEGA_DEV_CUT
        (30 deg, deviation of omega from nearest cis 0 / trans 180);
        mainchain_outlier when out-of-plane n_oop_deg >
        MAINCHAIN_OOP_CUT (10 deg). Former BACKBONE_GATES tolerates
        peptide_outlier_frac=0.05 / mainchain_outlier_frac=0.05.
        omega_dev vs nearest cis(0)/trans(180), outlier when > 30 deg;
        n_oop_deg of C(i)->N(i+1) out of CA(i)/C(i)/O(i) plane, outlier
        when > 10 deg.
        '''
        res = self.standard_residues()
        rows = []
        for i, r in enumerate(res):
            if i >= len(res) - 1:
                break
            nxt = res[i + 1]
            if not (all(a in r for a in ('CA', 'C', 'O'))
                    and all(a in nxt for a in ('N', 'CA'))):
                continue
            omega = self.dihedral_angle(r['CA'].coord, r['C'].coord,
                                        nxt['N'].coord, nxt['CA'].coord)
            dev = min(abs(omega), abs(omega - 180.0), abs(omega + 180.0))
            CA, C, O, Nn = (r['CA'].coord, r['C'].coord,
                            r['O'].coord, nxt['N'].coord)
            norm = np.cross(CA - C, O - C)
            norm /= np.linalg.norm(norm)
            cn = Nn - C
            cn /= np.linalg.norm(cn)
            oop = float(np.degrees(np.arcsin(
                np.clip(abs(np.dot(norm, cn)), 0, 1))))
            rows.append([r.id[1], r.resname, omega, dev,
                         dev > self.OMEGA_DEV_CUT,
                         oop, oop > self.MAINCHAIN_OOP_CUT])
        cols = ['res_pos', 'resname', 'omega', 'omega_dev',
                'peptide_outlier', 'n_oop_deg', 'mainchain_outlier']
        return pd.DataFrame(rows, columns=cols)

    # ------------------------------------------------------------------
    # residue-type aware Ramachandran (MolProbity style)
    # ------------------------------------------------------------------

    def phi_psi_series(self) -> list:
        '''per-residue phi/psi in degrees with neighbour context.

        Returns list of dicts (res_pos, resname, phi, psi, prev_resname);
        phi/psi are None at chain ends or when backbone atoms are missing.
        '''
        res = self.standard_residues()
        out = []
        for i, r in enumerate(res):
            prev = res[i - 1] if i > 0 else None
            nxt = res[i + 1] if i < len(res) - 1 else None
            phi = psi = None
            if prev is not None and all(a in r for a in ('N', 'CA', 'C')) \
                    and 'C' in prev:
                phi = self.dihedral_angle(prev['C'].coord, r['N'].coord,
                                          r['CA'].coord, r['C'].coord)
            if nxt is not None and all(a in r for a in ('N', 'CA', 'C')) \
                    and 'N' in nxt:
                psi = self.dihedral_angle(r['N'].coord, r['CA'].coord,
                                          r['C'].coord, nxt['N'].coord)
            prev_name = prev.resname if prev is not None else None
            out.append({'res_pos': r.id[1], 'resname': r.resname,
                        'phi': phi, 'psi': psi, 'prev_resname': prev_name})
        return out

    @staticmethod
    def rama_class(phi, psi, resname, prev_resname=None):
        '''Residue-type aware phi/psi label: favoured/allowed/outlier/NA.

        Cutoffs: analytic phi/psi boxes in degrees (Gly/Pro/pre-Pro/general
        sets below); no numeric confidence cut - use RamaGrid (favoured 98%
        / allowed 99.95% density cuts) for empirical grading.
        Special-cases Gly (no CB), Pro (phi near -60) and pre-Pro, which
        the single-box evaluate_ramachandran() ignores. Analytic boxes only.
        '''
        if phi is None or psi is None \
                or not (np.isfinite(phi) and np.isfinite(psi)):
            return 'NA'
        gly = resname == 'GLY'
        pro = resname == 'PRO'
        prepro = (prev_resname == 'PRO')

        def in_box(p0, p1, s0, s1):
            return p0 < phi < p1 and s0 <= psi < s1

        if gly:
            if in_box(-180, -30, 60, 180) or in_box(20, 180, -180, 0) \
                    or in_box(-180, 0, -90, 60) or in_box(0, 180, 60, 180):
                return 'favoured'
            return 'outlier' if in_box(0, 180, -180, -90) else 'allowed'
        if pro:
            if in_box(-90, -30, 60, 180) or in_box(-90, -30, -90, 60):
                return 'favoured'
            if not in_box(-110, -10, -180, 180):
                return 'outlier'
            return 'allowed'
        alpha = in_box(-110, -30, -75, 10)
        beta = in_box(-180, -40, 90, 180)
        lh = in_box(20, 90, -10, 60)
        if prepro:
            ext = in_box(-180, -40, 100, 180) \
                or in_box(-90, -20, -80, 100) or in_box(-90, 20, 20, 90)
            if ext or alpha or beta:
                return 'favoured'
            if in_box(-180, 40, -180, 180):
                return 'allowed'
            return 'outlier'
        if alpha or beta:
            return 'favoured'
        if in_box(-180, 0, -110, 180) or in_box(0, 120, -180, 60) \
                or lh or in_box(-120, -40, 10, 90):
            return 'allowed'
        return 'outlier'

    def rama_molprobity(self) -> dict:
        '''Residue-type aware Ramachandran summary from phi_psi_series().

        Cutoff: none on fractions (informational) - check() gates
        rama_outlier_frac (and torsion null fraction) only.
        Returns analysed count, favoured/allowed/outlier fractions and
        outliers list of (pos, resname, phi, psi).
        '''
        counts = {'favoured': 0, 'allowed': 0, 'outlier': 0, 'NA': 0}
        outliers = []
        for row in self.phi_psi_series():
            cls = self.rama_class(row['phi'], row['psi'], row['resname'],
                                  row['prev_resname'])
            counts[cls] += 1
            if cls == 'outlier':
                rec = (row['res_pos'], row['resname'], row['phi'], row['psi'])
                outliers.append(rec)
        analysed = sum(counts[k] for k in ('favoured', 'allowed', 'outlier'))
        frac = {k: (counts[k] / analysed if analysed else 0.0)
                for k in ('favoured', 'allowed', 'outlier')}
        return {
            'analysed': analysed,
            **{k: round(v, 3) for k, v in frac.items()},
            'outliers': outliers,
        }

    def rama_grid_summary(self, grids) -> dict:
        '''Empirical grid-based Ramachandran summary of one chain.

        grids: RamaGrid instance (structure_profiling/geometry/rama_grid.py).
        Cutoff: none on fractions (informational; same shape as
        Torsions.rama_molprobity()) - per-residue cuts are favoured 98% /
        allowed 99.95% inside RamaGrid.classify().
        '''
        counts = {'favoured': 0, 'allowed': 0, 'outlier': 0, 'NA': 0}
        outliers = []
        for row in self.phi_psi_series():
            label = grids.classify(row['phi'], row['psi'], row['resname'],
                                    row['prev_resname'])
            counts[label] += 1
            if label == 'outlier':
                rec = (row['res_pos'], row['resname'], row['phi'], row['psi'])
                outliers.append(rec)
        analysed = sum(counts[k] for k in ('favoured', 'allowed', 'outlier'))
        frac = {k: (counts[k] / analysed if analysed else 0.0)
                for k in ('favoured', 'allowed', 'outlier')}
        return {
            'analysed': analysed,
            **{k: round(v, 3) for k, v in frac.items()},
            'outliers': outliers,
        }

    # ------------------------------------------------------------------
    # summary + gates (used by pipelines/qc/chain_torsions.py)
    # ------------------------------------------------------------------

    def run(self, keep_details: bool = True) -> dict:
        '''compute all torsion / bond-angle metrics of one chain.

        Covers the torsion half of Geometry.run() and the planarity /
        Rama parts of the former structure_qc_screen():
        backbone bond angles vs Engh & Huber ideals (n_angles,
        angle_rmsz, n_z_outliers, z_outlier_frac), backbone torsions
        phi/psi/omega (n_torsions, torsion_null_frac, frac_trans,
        frac_cis, n_cis_nonpro), MolProbity-style Ramachandran
        (rama_analysed, rama_favoured_frac, rama_allowed_frac,
        rama_outlier_frac) and wwPDB peptide / main-chain planarity
        (n_peptides, peptide_outlier_frac, mainchain_outlier_frac).

        Omega fractions use the measurable peptides only (omega is None
        at fragment ends); torsion_null_frac is the fraction of residues
        whose phi/psi could not be measured (same quantity as the
        'torsion_angles' gate in Geometry.check()).

        Returns a dict of summary statistics matching the chain_torsions
        SQL table; raw tables are nested under 'details' when keep_details.
        '''
        res = self.standard_residues()
        angles = self.bond_angle_z()
        tors = self.torsion_angles()
        plan = self.peptide_planarity()
        rama = self.rama_molprobity()

        n_ang = len(angles)
        n_out = int(angles['outlier'].sum()) if n_ang else 0
        # omega statistics over the measurable peptides only
        om = tors[tors['omega'].notna()] if len(tors) else tors
        conf = om['conformation'].value_counts(normalize=True) \
            if len(om) else pd.Series(dtype=float)
        cis = om[om['conformation'] == 'cis'] if len(om) else om
        info = {
            'n_residues': len(res),
            # backbone bond angles (N-CA-C, CA-C-N, C-N-CA)
            'n_angles': n_ang,
            'angle_rmsz': round(
                float(np.sqrt(np.mean(angles['z'] ** 2))), 3
            ) if n_ang else 0.0,
            'n_z_outliers': n_out,
            'z_outlier_frac': round(n_out / n_ang, 3) if n_ang else 0.0,
            # backbone torsions phi / psi / omega
            'n_torsions': len(tors),
            'torsion_null_frac': round(
                float((tors['ramachandran'] == 'null').mean()), 3
            ) if len(tors) else 0.0,
            'frac_trans': round(float(conf.get('trans', 0.0)), 3),
            'frac_cis': round(float(conf.get('cis', 0.0)), 3),
            'n_cis_nonpro': int((cis['residue'] != 'PRO').sum())
                if len(cis) else 0,
            # residue-type aware Ramachandran (MolProbity style)
            'rama_analysed': rama['analysed'],
            'rama_favoured_frac': rama['favoured'],
            'rama_allowed_frac': rama['allowed'],
            'rama_outlier_frac': rama['outlier'],
            # wwPDB peptide / main-chain planarity
            'n_peptides': len(plan),
            'peptide_outlier_frac': round(
                float(plan['peptide_outlier'].sum()) / len(plan), 3
            ) if len(plan) else 0.0,
            'mainchain_outlier_frac': round(
                float(plan['mainchain_outlier'].sum()) / len(plan), 3
            ) if len(plan) else 0.0,
        }
        if keep_details:
            info['details'] = {
                'bond_angle_z': angles,
                'torsion_angles': tors,
                'peptide_planarity': plan,
                'rama_outliers': rama['outliers'],
            }
        return info

    def check(self, info: dict = None) -> tuple:
        '''gate on the summary statistics returned by run().

        Cutoffs (fail -> reason), see Torsions.GATES:
          angle_rmsz > 1.5          -> 'angle_rmsz' (RMS |Z| of angles);
          z_outlier_frac > 0.01     -> 'z_outliers' (|Z| > 5 angles);
          torsion_null_frac > 0.05  -> 'torsion_null' (same gate as
                                        Geometry.check 'torsion_angles');
          peptide_outlier_frac > 0.05 -> 'peptide_planarity'
                                        (omega dev > 30 deg);
          mainchain_outlier_frac > 0.05 -> 'mainchain_planarity'
                                        (out-of-plane > 10 deg);
          rama_outlier_frac > 0.05  -> 'rama_outliers' (MolProbity style).
        '''
        if info is None:
            info = self.run(keep_details=False)
        if not info:
            return False, 'null_torsions'
        gates = self.GATES
        # backbone bond angles: RMS-Z magnitude and |Z| > 5 outlier fraction
        if info.get('angle_rmsz', 0.0) > gates['angle_rmsz']:
            return False, 'angle_rmsz'
        if info.get('z_outlier_frac', 0.0) > gates['z_outlier_frac']:
            return False, 'z_outliers'
        # residues without measurable phi/psi (same gate as Geometry.check)
        if info.get('torsion_null_frac', 0.0) > gates['torsion_null_frac']:
            return False, 'torsion_null'
        # wwPDB planarity: omega dev > 30 deg / out-of-plane > 10 deg
        if info.get('peptide_outlier_frac', 0.0) > \
                gates['peptide_outlier_frac']:
            return False, 'peptide_planarity'
        if info.get('mainchain_outlier_frac', 0.0) > \
                gates['mainchain_outlier_frac']:
            return False, 'mainchain_planarity'
        # MolProbity-style Ramachandran outliers
        if info.get('rama_outlier_frac', 0.0) > gates['rama_outlier_frac']:
            return False, 'rama_outliers'
        return True, 'pass'