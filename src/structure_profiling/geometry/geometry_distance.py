'''Distance-based protein geometry QC for one Bio.PDB chain.

Covers standard residues / heavy atoms, radius of gyration, covalent bond
length Z-scores (Engh & Huber style, |Z| > 5 outlier) and MolProbity-style
Cb deviation.

Clash analysis (steric clashes, MolProbity-style vdW-overlap clash score)
is collected in structure/clashes.py (Clashes); the methods and
module-level functions below delegate to it for backward compatibility.
Chain-continuity analysis (missing backbone atoms, chain breaks,
Ca-Ca / peptide C-N distances) is collected in qc/breaks.py (Breaks);
the methods below delegate to it for backward compatibility.
'''
import pandas as pd
from Bio import PDB
import numpy as np
from pathlib import Path

from ..qc import Clashes

class GeometryDistance:
    # canonical clash cutoffs live in Clashes (structure/clashes.py)
    severe_clash_distance:float = Clashes.severe_clash_distance

    # ideal backbone bond lengths (mean, sd) in A, Engh & Huber style (cf.
    # wwPDB standard-geometry compilation). Keys are sorted atom-name pairs.
    IDEAL_BONDS = {
        ('CA', 'N'): (1.458, 0.019),
        ('C', 'CA'): (1.525, 0.021),
        ('C', 'N'): (1.329, 0.014),    # peptide bond
        ('C', 'O'): (1.231, 0.020),
        ('CA', 'CB'): (1.521, 0.033),
    }
    CB_BOND_LEN = 1.521      # CA-CB (A)
    CB_ANGLE = 110.0         # N-CA-CB = C-CA-CB (deg), idealized
    CB_DEVIATION_CUT = 0.25  # A, MolProbity-style Cb outlier threshold
    Z_OUTLIER_CUT = 5.0      # |Z| above this flags a covalent outlier

    def __init__(self, chain):
        self.chain = chain

    def standard_residues(self):
        '''standard amino-acid residues (skip hetero / water).

        Cutoff: none - filters on residue id[0] == ' ' only.
        '''
        return [r for r in self.chain if r.id[0] == " "]

    def heavy_atoms(self):
        '''heavy (non-hydrogen) atoms of standard residues.

        Cutoff: none - filters on element != 'H' only.
        '''
        return [
            a for r in self.standard_residues()
            for a in r if a.element != 'H'
        ]

    def missing_residues(self) -> list:
        '''check if every residues has expected backbone atoms: N, CA, C, O
        return: residues with missing atoms

        Moved to Breaks.missing_residues (qc/breaks.py); this delegate
        keeps the GeometryDistance API used by Geometry.run().

        Cutoff: no tolerance - any residue lacking N/CA/C/O is listed;
        Geometry.check() / Breaks.check() fail when num_missing > 0.
        '''
        from ..qc.breaks import Breaks
        return Breaks(self.chain).missing_residues()

    def chain_continuity(self, cn_break: float = 2.0,
            ca_break: float = 4.5, cn_strain: float = 1.5):
        '''Localise backbone breaks / chain-continuity failures for one chain.

        Moved to Breaks.chain_continuity (qc/breaks.py); this delegate
        keeps the GeometryDistance API used by Geometry.run().

        Cutoff: cn_break=2.0 A, C-N above this is a real break (ideal
        1.329 A); ca_break=4.5 A, CA-CA above this is a break (ideal
        ~3.8 A); cn_strain=1.5 A, C-N above this but <= cn_break is a
        stretched peptide. A pair is also a break when C-N is unmeasurable
        (missing atoms) or residue numbering skips (resnum_gap != 1).
        '''
        from ..qc.breaks import Breaks
        return Breaks(self.chain).chain_continuity(
            cn_break=cn_break, ca_break=ca_break, cn_strain=cn_strain)

    def ca_distances(self):
        '''detect bad Ca connections. No large gap
        calculate distance between Ca-Ca atoms

        Moved to Breaks.ca_distances (qc/breaks.py); this delegate keeps
        the GeometryDistance API used by Geometry.run().

        Cutoff: ideal CA-CA ~3.8 A; Geometry.check() fails when max > 5 A.
        '''
        from ..qc.breaks import Breaks
        return Breaks(self.chain).ca_distances()

    def peptide_distances(self):
        """peptide C-N distance
        https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111178&ver=0&trm=peptide&org=

        Moved to Breaks.peptide_distances (qc/breaks.py); this delegate
        keeps the GeometryDistance API used by Geometry.run().

        Cutoff: ideal C-N = 1.329 A (peptide bond); chain_continuity()
        flags C-N > 2.0 A as a break; Geometry.check() fails when
        max > 5 A or mean > 3 A.
        """
        from ..qc.breaks import Breaks
        return Breaks(self.chain).peptide_distances()

    def atom_distance(self, pos1:int, atom1:str, pos2:int, atom2:str):
        '''
        distance (A) between two named atoms of one chain, e.g.
        atom_distance(chain, 40, 'SG', 114, 'SG') for a disulfide.
        returns None if chain / residue / atom is missing.
        Cutoff: none - raw distance in A (disulfide SG-SG ideal ~2.05 A).
        '''
        residues = {r.id[1]: r for r in self.standard_residues()}
        if pos1 not in residues or pos2 not in residues:
            return None
        r1, r2 = residues[pos1], residues[pos2]
        if atom1 not in r1 or atom2 not in r2:
            return None
        return float(np.linalg.norm(r1[atom1].coord - r2[atom2].coord))

    def radius_of_gyration(self, all_atoms:bool=False):
        '''compactness: radius of gyration (A) of Ca atoms, or all heavy atoms.

        Cutoff: none (informational only) - reported rounded to 2 dp in
        Geometry.run(); not gated by check().
        '''
        if all_atoms:
            coords = [a.coord for a in self.heavy_atoms()]
        else:
            coords = [r["CA"].coord for r in self.standard_residues() if "CA" in r]
        ca = np.array(coords, dtype=float)
        if len(ca) < 2:
            return None
        center = ca.mean(axis=0)
        return float(np.sqrt(((ca - center) ** 2).sum(axis=1).mean()))

    def steric_clashes(self, cutoff:float=2.4):
        '''
        steric clashes: non-bonded heavy-atom contacts closer than cutoff (A).

        Moved to Clashes.steric_clashes (structure/clashes.py); this
        delegate keeps the GeometryDistance API used by Geometry.run().

        Cutoffs: cutoff=2.4 A fixed distance (heuristic); contacts below
        severe_clash_distance=2.0 A are counted severe (Clashes.severe_clashes);
        clash_score (per 1000 atoms) fails Clashes.check() above 30.0.

        covalently bonded atoms are not clashes: pairs inside one residue,
        between sequential residues (i, i+1) and disulfide SG-SG bonds are
        excluded. hydrogens and hetero / water atoms are ignored.
        returns DataFrame(atom_pos, atom, other_pos, other_atom, distance).
        '''
        return Clashes(self.chain).steric_clashes(cutoff)

    @staticmethod
    def _bond_z(name1:str, coord1, name2:str, coord2):
        '''Z-score of the |coord1 - coord2| distance against IDEAL_BONDS.

        Cutoff: none here - caller flags |Z| > Z_OUTLIER_CUT (5.0).
        Returns (z, observed, ideal, sd), or None for unknown atom pairs.
        '''
        key = tuple(sorted((name1, name2)))
        if key not in GeometryDistance.IDEAL_BONDS:
            return None
        mean, sd = GeometryDistance.IDEAL_BONDS[key]
        obs = float(np.linalg.norm(coord1 - coord2))
        return (obs - mean) / sd, obs, mean, sd

    def bond_length_z(self) -> pd.DataFrame:
        '''
        Z-scores of backbone bond lengths against IDEAL_BONDS (|Z| > 5 is
        an outlier, wwPDB criterion).

        Cutoff: Z_OUTLIER_CUT=5.0 on |Z|; bond RMS-Z expected in [0, 1]
        (wwPDB), former BACKBONE_GATES tolerates bond_rmsz up to 2.0.

        Evaluates intra-residue N-CA, CA-C, C-O, CA-CB plus the C(i)-N(i+1)
        peptide bridge. Returns DataFrame(res_pos, resname, atoms, observed,
        ideal, sd, z, outlier).
        '''
        cols = ['res_pos', 'resname', 'atoms', 'observed',
                'ideal', 'sd', 'z', 'outlier']
        rows = []
        res = self.standard_residues()
        for i, r in enumerate(res):
            for key in (('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB')):
                if all(a in r for a in key):
                    scored = self._bond_z(key[0], r[key[0]].coord,
                                          key[1], r[key[1]].coord)
                    if scored is not None:
                        z, obs, mean, sd = scored
                        rows.append([r.id[1], r.resname, '-'.join(key),
                                     obs, mean, sd, z,
                                     abs(z) > self.Z_OUTLIER_CUT])
            if i < len(res) - 1 and 'C' in r and 'N' in res[i + 1]:
                scored = self._bond_z('C', r['C'].coord,
                                      'N', res[i + 1]['N'].coord)
                if scored is not None:
                    z, obs, mean, sd = scored
                    rows.append([r.id[1], r.resname, 'C-N',
                                 obs, mean, sd, z,
                                 abs(z) > self.Z_OUTLIER_CUT])
        return pd.DataFrame(rows, columns=cols)

    @staticmethod
    def virtual_cb(N, CA, C):
        '''
        Idealized Cb position for an L-amino acid built from the N, CA, C
        frame (Lovell et al. 2003 style).

        Cutoffs: CA-CB = CB_BOND_LEN (1.521 A); N-CA-CB = C-CA-CB =
        CB_ANGLE (110.0 deg); cbeta_deviation() flags dev > 0.25 A.
        '''
        e1 = (N - CA) / np.linalg.norm(N - CA)
        e2 = (C - CA) / np.linalg.norm(C - CA)
        n = np.cross(e1, e2)
        n /= np.linalg.norm(n)
        cosA = np.cos(np.radians(GeometryDistance.CB_ANGLE))
        d = float(np.dot(e1, e2))
        a, b = np.linalg.solve([[1.0, d], [d, 1.0]], [cosA, cosA])
        c2 = 1.0 - a * a - b * b - 2.0 * a * b * d
        c = np.sqrt(max(c2, 0.0))
        u = a * e1 + b * e2 + c * n
        return CA + GeometryDistance.CB_BOND_LEN * u

    def cbeta_deviation(self) -> pd.DataFrame:
        '''
        MolProbity-style Cb deviation (A): distance between the observed CB
        and the idealized position from virtual_cb().

        Cutoff: CB_DEVIATION_CUT=0.25 A on cb_dev (MolProbity-style).

        Returns DataFrame(res_pos, resname, cb_dev, outlier
        with cb_dev > CB_DEVIATION_CUT). Gly residues (no CB) report
        cb_dev None and outlier False. Large deviations mark backbone
        errors - small on clean models, huge on chain-broken ones.
        '''
        rows = []
        for r in self.standard_residues():
            if not all(a in r for a in ('N', 'CA', 'C')) or 'CB' not in r:
                rows.append([r.id[1], r.resname, None, False])
                continue
            pred = self.virtual_cb(r['N'].coord, r['CA'].coord, r['C'].coord)
            dev = float(np.linalg.norm(pred - r['CB'].coord))
            rows.append([r.id[1], r.resname, dev, dev > self.CB_DEVIATION_CUT])
        return pd.DataFrame(rows, columns=['res_pos', 'resname', 'cb_dev',
                                           'outlier'])

    def vdw_clashes(self, radii=None, overlap_cut: float = 0.40,
                    include_hydrogen: bool = True) -> pd.DataFrame:
        '''MolProbity-style vdW-overlap clashes for this chain.

        Moved to Clashes.vdw_clashes (structure/clashes.py); delegate.

        Cutoff: overlap_cut=0.40 A on vdW overlap (r_i + r_j - d_ij;
        VDW_OVERLAP_CUT, MolProbity definition); radii from VDW_RADII.
        '''
        return Clashes(self.chain).vdw_clashes(radii=radii,
                                               overlap_cut=overlap_cut,
                                               include_hydrogen=include_hydrogen)

    def vdw_clash_score(self, clashes: pd.DataFrame,
                        n_atoms: int = None) -> dict:
        '''Clash score (per 1000 atoms) for a vdw_clashes() table.

        Moved to Clashes.vdw_clash_score (structure/clashes.py); delegate.

        Cutoff: MolProbity clashscore bands (<5 very good, 5-10 acceptable,
        10-20 investigate, >20 poor); Geometry.check() analogue fails above
        clash_score_cutoff=30.0.
        '''
        return Clashes(self.chain).vdw_clash_score(clashes, n_atoms=n_atoms)


# ===========================================================================
# MolProbity-style van-der-Waals clash score
# (implemented in structure/clashes.py; aliases kept for backward compat)
# ===========================================================================

VDW_RADII = Clashes.VDW_RADII
VDW_OVERLAP_CUT = Clashes.VDW_OVERLAP_CUT


def add_hydrogens_pdbfixer(pdb_file, out_file=None, ph=7.0):
    '''Add missing heavy atoms + hydrogens with PDBFixer/OpenMM.'''
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=str(pdb_file))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]      # ignore terminal gaps
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    out_file = out_file or str(Path(pdb_file).with_suffix('_h.pdb'))
    with open(out_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    return Path(out_file)


def vdw_clashes(chain, radii=None, overlap_cut: float = VDW_OVERLAP_CUT,
                include_hydrogen: bool = True) -> pd.DataFrame:
    '''vdW-overlap clashes for one chain (MolProbity definition).

    Moved to Clashes.vdw_clashes (structure/clashes.py); module-level
    delegate keeps the old function API.

    Cutoff: overlap_cut=0.40 A (VDW_OVERLAP_CUT) on overlap =
    (r_i + r_j) - d_ij; same-residue, sequential (i, i+1) and SG-SG
    disulfide pairs excluded.
    '''
    return Clashes(chain).vdw_clashes(radii=radii, overlap_cut=overlap_cut,
                                      include_hydrogen=include_hydrogen)


def vdw_clash_score(clashes: pd.DataFrame, chain, n_atoms: int = None):
    '''clashes per 1000 atoms and overlap stats (MolProbity clashscore).

    Moved to Clashes.vdw_clash_score (structure/clashes.py); delegate.

    Cutoff: clashscore bands <5 very good / 5-10 acceptable / 10-20
    investigate / >20 poor; 0 means no serious overlaps found.
    '''
    return Clashes(chain).vdw_clash_score(clashes, n_atoms=n_atoms)
