'''Angle-based protein geometry QC for one Bio.PDB chain.

Covers bond angles, dihedral angles, single-box Ramachandran evaluation,
bond-angle Z-scores (Engh & Huber style), CA chirality (L/D), peptide and
main-chain planarity (wwPDB criteria) and residue-type aware Ramachandran
classification (MolProbity style: Gly / Pro / pre-Pro special cases).
Empirical Top8000-style grids live in geometry_rama.py.
'''
import numpy as np
import pandas as pd
from Bio import PDB

from .rama_grid import RamaGrid

class GeometryAngle:

    # ideal backbone bond angles (mean, sd) in degrees, Engh & Huber style
    IDEAL_ANGLES = {
        ('N', 'CA', 'C'): (111.2, 2.8),
        ('CA', 'C', 'N'): (116.2, 2.0),
        ('C', 'N', 'CA'): (121.9, 2.6),
    }
    CHIRAL_VOL_L = 2.6    # A^3, L-amino-acid CA chiral volume
    CHIRAL_VOL_MIN = 0.7  # below this the CA centre is 'planar' / broken
    Z_OUTLIER_CUT = 5.0   # |Z| above this flags a covalent outlier

    def __init__(self, chain):
        self.chain = chain

    @staticmethod
    def dihedral_angle(p0, p1, p2, p3):
        '''dihedral angle (degrees) between four points p0-p1-p2-p3'''
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

    @staticmethod
    def bond_angle(a, b, c):
        '''bond angle (degrees) a-b-c around the central atom b'''
        v1, v2 = a - b, c - b
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
    
    @staticmethod
    def evaluate_ramachandran(phi, psi):
        """
        Check if phi/psi angles are in allowed regions
        General allowed regions for L-amino acids (degrees)
        https://jenalib.leibniz-fli.de/ImgLibDoc/glossary/IMAGE_RAMA.htm

        Cutoffs (analytic boxes, degrees): beta-sheet phi in (-180, -60)
        and psi in (90, 180); right-handed helix phi in (-120, -30), psi
        in (-60, 30); left-handed helix phi in (45, 90), psi in (0, 60);
        'allowed' phi in [-180, 0], psi in [0, 180]; else 'bad'.
        phi/psi are expected in degrees; None / NaN -> 'null'
        """
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

    def torsion_angles(self):
        '''backbone torsion angles (degrees)
        Calculate: φ (phi), ψ (psi), ω (omega) and look for Ramachandran
        outliers. omega is None at the end of each peptide fragment;
        conformation is cis (|omega| < 90 deg) or trans.
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
                        conformation = 'cis' if abs(omega) < 90 else 'trans'
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

    @staticmethod
    def _std(chain):
        '''standard amino-acid residues (skip hetero / water)'''
        return [r for r in chain if r.id[0] == ' ']

    def bond_angle_z(self) -> pd.DataFrame:
        '''Z-scores of backbone bond angles vs IDEAL_ANGLES (|Z| > 5 outlier).

        Cutoff: Z_OUTLIER_CUT=5.0 on |Z|; angle RMS-Z expected in [0, 1]
        (wwPDB), BACKBONE_GATES tolerates angle_rmsz up to 1.5 for design.
        Covers intra-residue N-CA-C plus CA(i)-C(i)-N(i+1) and
        C(i)-N(i+1)-CA(i+1) bridge angles. Returns columns res_pos, resname,
        atoms, observed, ideal, sd, z, outlier.
        '''
        cols = ['res_pos', 'resname', 'atoms', 'observed',
                'ideal', 'sd', 'z', 'outlier']
        rows = []
        res = self._std(self.chain)
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

    def ca_chirality(self) -> pd.DataFrame:
        '''CA chiral volume from N, CA, C, CB (~ +2.6 A^3 L, -2.6 D).

        Cutoffs: CHIRAL_VOL_L=2.6 A^3 expected L volume; |vol| below
        CHIRAL_VOL_MIN=0.7 A^3 is 'planar' (broken); sign vs majority
        reference decides L/D.
        Returns res_pos, resname, chiral_volume, handedness (L/D/planar/NA).
        Reference sign is the majority; |vol| < CHIRAL_VOL_MIN is planar;
        residues without N/CA/C/CB (e.g. Gly) are NA.
        '''
        rows = []
        for r in self._std(self.chain):
            if not all(a in r for a in ('N', 'CA', 'C', 'CB')):
                rows.append([r.id[1], r.resname, None, 'NA'])
                continue
            N, CA, C, CB = (r[a].coord for a in ('N', 'CA', 'C', 'CB'))
            vol = float(np.dot(CB - CA, np.cross(N - CA, C - CA)))
            rows.append([r.id[1], r.resname, vol, None])
        df = pd.DataFrame(rows, columns=['res_pos', 'resname',
                                         'chiral_volume', 'handedness'])
        sign = np.sign(df['chiral_volume'].dropna())
        if len(sign) == 0:
            return df
        ref = sign.median()
        df['handedness'] = np.where(
            df['chiral_volume'].isna(), 'NA',
            np.where(df['chiral_volume'].abs() < self.CHIRAL_VOL_MIN, 'planar',
                     np.where(np.sign(df['chiral_volume']) == ref,
                              'L' if ref > 0 else 'D',
                              'D' if ref > 0 else 'L')))
        return df

    def peptide_planarity(self) -> pd.DataFrame:
        '''wwPDB planarity per peptide unit (i -> i+1).

        Cutoffs (wwPDB): peptide_outlier when omega_dev > 30 deg (deviation
        of omega from nearest cis 0 / trans 180); mainchain_outlier when
        out-of-plane n_oop_deg > 10 deg. BACKBONE_GATES tolerates
        peptide_outlier_frac=0.05 / mainchain_outlier_frac=0.05.
        omega_dev vs nearest cis(0)/trans(180), outlier when > 30 deg;
        n_oop_deg of C(i)->N(i+1) out of CA(i)/C(i)/O(i) plane, outlier
        when > 10 deg.
        '''
        res = self._std(self.chain)
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
            rows.append([r.id[1], r.resname, omega, dev, dev > 30,
                         oop, oop > 10])
        cols = ['res_pos', 'resname', 'omega', 'omega_dev',
                'peptide_outlier', 'n_oop_deg', 'mainchain_outlier']
        return pd.DataFrame(rows, columns=cols)

    def phi_psi_series(self) -> list:
        '''per-residue phi/psi in degrees with neighbour context.

        Returns list of dicts (res_pos, resname, phi, psi, prev_resname);
        phi/psi are None at chain ends or when backbone atoms are missing.
        '''
        res = self._std(self.chain)
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
            return p0 <= phi < p1 and s0 <= psi < s1

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

        Cutoff: none on fractions (informational) - Geometry.check()
        analogue gates torsion null fraction > 0.05 only.
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

    def rama_grid_summary(self, grids: RamaGrid) -> dict:
        '''Empirical grid-based Ramachandran summary of one chain.

        Cutoff: none on fractions (informational; same shape as
        GeometryAngle.rama_molprobity()) - per-residue cuts are favoured 98% /
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
