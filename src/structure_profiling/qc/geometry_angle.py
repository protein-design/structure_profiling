'''Angle-based protein geometry QC for one Bio.PDB chain.

Covers CA chirality (L/D) directly; bond angles, dihedral angles,
single-box Ramachandran evaluation, bond-angle Z-scores (Engh & Huber
style), peptide / main-chain planarity (wwPDB criteria) and residue-type
aware Ramachandran classification (MolProbity style: Gly / Pro / pre-Pro)
are collected in qc/torsions.py (Torsions); the methods below delegate
to it for backward compatibility (lazy imports - torsions.py imports
this module at load time). Empirical Top8000-style grids live in
rama_grid.py.

CA chiral-volume cutoffs (from wwPDB validation):
    CHIRAL_VOL_L = 2.6    ideal L-amino-acid CA chiral volume (A^3)
    CHIRAL_VOL_MIN = 0.7  below this the CA centre is 'planar'/broken
The angle cutoffs (IDEAL_ANGLES, Z_OUTLIER_CUT) are defined here and
re-exported by Torsions (qc/torsions.py).
'''
import numpy as np
import pandas as pd


class GeometryAngle:

    # ideal backbone bond angles (mean, sd) in degrees, Engh & Huber style;
    # canonical values re-exported by Torsions (qc/torsions.py)
    IDEAL_ANGLES = {
        ('N', 'CA', 'C'): (111.2, 2.8),
        ('CA', 'C', 'N'): (116.2, 2.0),
        ('C', 'N', 'CA'): (121.9, 2.6),
    }
    Z_OUTLIER_CUT: float = 5.0  # |Z| above this flags a covalent outlier
    CHIRAL_VOL_L = 2.6    # A^3, L-amino-acid CA chiral volume
    CHIRAL_VOL_MIN = 0.7  # below this the CA centre is 'planar' / broken

    def __init__(self, chain):
        self.chain = chain

    @staticmethod
    def _std(chain):
        '''standard amino-acid residues (skip hetero / water)'''
        return [r for r in chain if r.id[0] == ' ']

    # ------------------------------------------------------------------
    # delegates to Torsions (qc/torsions.py)
    # ------------------------------------------------------------------

    @staticmethod
    def dihedral_angle(p0, p1, p2, p3):
        '''dihedral angle (degrees) between four points p0-p1-p2-p3.

        Moved to Torsions.dihedral_angle (qc/torsions.py); delegate keeps
        the GeometryAngle API used by geometry.py and the tests.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions.dihedral_angle(p0, p1, p2, p3)

    @staticmethod
    def bond_angle(a, b, c):
        '''bond angle (degrees) a-b-c around the central atom b.

        Moved to Torsions.bond_angle (qc/torsions.py); delegate.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions.bond_angle(a, b, c)

    @staticmethod
    def evaluate_ramachandran(phi, psi):
        '''
        Check if phi/psi angles are in allowed regions
        General allowed regions for L-amino acids (degrees)
        https://jenalib.leibniz-fli.de/ImgLibDoc/glossary/IMAGE_RAMA.htm

        Moved to Torsions.evaluate_ramachandran (qc/torsions.py);
        delegate keeps the API used by geometry.py and the tests.
        phi/psi are expected in degrees; None / NaN -> 'null'
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions.evaluate_ramachandran(phi, psi)

    def bond_angles(self):
        '''raw backbone bond angles (N-CA-C, CA-C-N; ideal ~111 / ~117 deg).

        Moved to Torsions.bond_angles (qc/torsions.py); delegate.
        Returns DataFrame(residue, res_pos, angle_N_CA_C, angle_CA_C_N);
        angle_CA_C_N is None at peptide fragment ends.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).bond_angles()

    def torsion_angles(self):
        '''backbone torsion angles (degrees): phi / psi / omega per residue
        with cis/trans conformation (|omega| < 90 deg is cis) and
        single-box Ramachandran label.

        Moved to Torsions.torsion_angles (qc/torsions.py); delegate.
        omega is None at the end of each peptide fragment.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).torsion_angles()

    def bond_angle_z(self) -> pd.DataFrame:
        '''Z-scores of backbone bond angles vs IDEAL_ANGLES (|Z| > 5 outlier).

        Moved to Torsions.bond_angle_z (qc/torsions.py); delegate.
        Covers intra-residue N-CA-C plus CA(i)-C(i)-N(i+1) and
        C(i)-N(i+1)-CA(i+1) bridge angles. Returns columns res_pos,
        resname, atoms, observed, ideal, sd, z, outlier.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).bond_angle_z()

    def peptide_planarity(self) -> pd.DataFrame:
        '''wwPDB planarity per peptide unit (i -> i+1).

        Moved to Torsions.peptide_planarity (qc/torsions.py); delegate.
        omega_dev vs nearest cis(0)/trans(180), outlier when > 30 deg;
        n_oop_deg of C(i)->N(i+1) out of CA(i)/C(i)/O(i) plane, outlier
        when > 10 deg.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).peptide_planarity()

    def phi_psi_series(self) -> list:
        '''per-residue phi/psi in degrees with neighbour context.

        Moved to Torsions.phi_psi_series (qc/torsions.py); delegate.
        Returns list of dicts (res_pos, resname, phi, psi, prev_resname);
        phi/psi are None at chain ends or when backbone atoms are missing.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).phi_psi_series()

    @staticmethod
    def rama_class(phi, psi, resname, prev_resname=None):
        '''Residue-type aware phi/psi label: favoured/allowed/outlier/NA.

        Moved to Torsions.rama_class (qc/torsions.py); delegate keeps the
        API used by the tests. Analytic boxes only.
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions.rama_class(phi, psi, resname, prev_resname)

    def rama_molprobity(self) -> dict:
        '''Residue-type aware Ramachandran summary from phi_psi_series().

        Moved to Torsions.rama_molprobity (qc/torsions.py); delegate.
        Returns analysed count, favoured/allowed/outlier fractions and
        outliers list of (pos, resname, phi, psi).
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).rama_molprobity()

    def rama_grid_summary(self, grids) -> dict:
        '''Empirical grid-based Ramachandran summary of one chain.

        Moved to Torsions.rama_grid_summary (qc/torsions.py); delegate.
        grids: RamaGrid instance (structure_profiling/qc/rama_grid.py).
        '''
        from .torsions import Torsions  # lazy, avoids import cycle
        return Torsions(self.chain).rama_grid_summary(grids)

    # ------------------------------------------------------------------
    # CA chirality (kept here: side-chain chirality, not a torsion/angle)
    # ------------------------------------------------------------------

    def ca_chirality(self) -> pd.DataFrame:
        '''CA chirality per standard residue via the determinant volume of
        N-CA-C-N signed by CA-CB (chiral volume, A^3).

        Cutoffs (wwPDB): volume > 0 is L (the CHIRAL_VOL_L=2.6 A^3 ideal
        for L-amino acids); 0 < volume < CHIRAL_VOL_MIN=0.7 A^3 flags a
        'planar' CA centre; volume < 0 is a D-amino acid (rare in natural
        proteins); residues without CB (glycine) are 'NA'. Returns
        columns res_pos, resname, chiral_vol, handedness.
        '''
        rows = []
        for r in self._std(self.chain):
            if 'CB' not in r or not all(a in r for a in ('N', 'CA', 'C')):
                rows.append([r.id[1], r.resname, None, 'NA'])
                continue
            N, CA, C, CB = (r['N'].coord, r['CA'].coord,
                            r['C'].coord, r['CB'].coord)
            vol = np.dot(np.cross(N - CA, C - CA), CB - CA)
            if vol < 0:
                chir = 'D'
            elif vol < self.CHIRAL_VOL_MIN:
                chir = 'planar'
            else:
                chir = 'L'
            rows.append([r.id[1], r.resname, float(vol), chir])
        cols = ['res_pos', 'resname', 'chiral_vol', 'handedness']
        return pd.DataFrame(rows, columns=cols)
