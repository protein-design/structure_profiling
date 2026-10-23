'''
Steric and van-der-Waals clash analysis for one Bio.PDB chain.

Collected here from src/structure_profiling (previously spread over
geometry/geometry_distance.py, geometry/geometry.py and qc/qc_screen.py):
- steric_clashes:   fixed-distance non-bonded heavy-atom contacts
- severe_clashes:   subset below severe_clash_distance (2.0 A)
- clash_score:      clashes per 1000 heavy atoms (heuristic)
- vdw_clashes:      MolProbity-style vdW overlaps (radii aware,
                    overlap >= VDW_OVERLAP_CUT = 0.40 A)
- vdw_clash_score:  vdW overlaps per 1000 heavy atoms (MolProbity clashscore)

run() returns the summary exported by pipelines/qc/chain_clashes.py to
the SQL table chain_clashes; check() gates it on the Clashes cutoffs.

references:
MolProbity (Chen et al. 2010, Read et al. 2011) and wwPDB OneDep
validation reports (X-ray + EM guides, updated 2025).
'''
import numpy as np
import pandas as pd
from Bio import PDB


class Clashes:
    '''clash analysis of one Bio.PDB chain (steric + MolProbity vdW).'''

    # steric clash QC thresholds
    clash_cutoff:float = 2.4           # A; contacts below are clashes
    severe_clash_distance:float = 2.0  # A; contacts below are severe clashes
    clash_score_cutoff:float = 30.0    # max clashes per 1000 heavy atoms
    # MolProbity clashscore bands (vdw_clash_score): <5 very good,
    # 5-10 acceptable, 10-20 investigate, >20 poor
    vdw_clash_cutoff:float = 20.0

    # non-hydrogen vdW radii (A); H (1.0 A) scored only when present / added
    VDW_RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
                 'P': 1.80, 'SE': 1.90, 'H': 1.00}
    VDW_OVERLAP_CUT = 0.40      # A; MolProbity flags overlap >= 0.4 A

    def __init__(self, chain):
        self.chain = chain

    # ------------------------------------------------------------------
    # atom selection (same filters as GeometryDistance)
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # steric clashes: fixed-distance non-bonded heavy-atom contacts
    # ------------------------------------------------------------------

    def steric_clashes(self, cutoff:float=None) -> pd.DataFrame:
        '''steric clashes: non-bonded heavy-atom contacts closer than
        cutoff (A).

        Cutoffs: cutoff=2.4 A (clash_cutoff) fixed distance (heuristic);
        contacts below severe_clash_distance=2.0 A are counted severe in
        severe_clashes(); clash_score (per 1000 atoms) fails check()
        above clash_score_cutoff=30.0.

        covalently bonded atoms are not clashes: pairs inside one residue,
        between sequential residues (i, i+1) and disulfide SG-SG bonds are
        excluded. hydrogens and hetero / water atoms are ignored.
        returns DataFrame(atom_pos, atom, other_pos, other_atom, distance).
        '''
        cutoff = self.clash_cutoff if cutoff is None else cutoff
        columns = ['atom_pos', 'atom', 'other_pos', 'other_atom', 'distance']
        atoms = self.heavy_atoms()
        if len(atoms) < 2:
            return pd.DataFrame(columns=columns)
        ns = PDB.NeighborSearch(atoms)
        seen, clashes = set(), []
        for atom in atoms:
            for other in ns.search(atom.coord, cutoff, level="A"):
                if atom is other:
                    continue
                key = tuple(sorted((
                    atom.get_serial_number(), other.get_serial_number()
                )))
                if key in seen:
                    continue
                seen.add(key)
                pos = atom.get_full_id()[3][1]
                other_pos = other.get_full_id()[3][1]
                if pos == other_pos:
                    continue                      # 1-2 / 1-3 bonds, same residue
                if abs(pos - other_pos) == 1:
                    continue                      # peptide-bond neighbours
                if atom.get_name() == 'SG' and other.get_name() == 'SG':
                    continue                      # disulfide bond
                clashes.append({
                    'atom_pos': pos,
                    'atom': atom.get_full_id()[4][0],
                    'other_pos': other_pos,
                    'other_atom': other.get_full_id()[4][0],
                    'distance': float(np.linalg.norm(atom.coord - other.coord)),
                })
        return pd.DataFrame(clashes, columns=columns)

    def severe_clashes(self, clashes:pd.DataFrame=None,
                       distance:float=None) -> pd.DataFrame:
        '''severe subset of a steric_clashes() table: contacts below
        distance (A).

        Cutoff: distance=severe_clash_distance (2.0 A); check() fails when
        num_severe_clashes > 0 (former BACKBONE_GATES num_servere_clashes).
        '''
        clashes = self.steric_clashes() if clashes is None else clashes
        distance = self.severe_clash_distance if distance is None else distance
        if clashes is None or clashes.empty:
            return clashes
        return clashes[clashes['distance'] < distance]

    def clash_score(self, clashes:pd.DataFrame=None, n_atoms:int=None) -> float:
        '''steric clashes per 1000 heavy atoms (heuristic score).

        Cutoff: none here - check() fails above clash_score_cutoff=30.0;
        reported rounded to 2 dp (consistent with vdw_clash_score()).
        '''
        clashes = self.steric_clashes() if clashes is None else clashes
        heavy = len(self.heavy_atoms()) if n_atoms is None else n_atoms
        if heavy == 0:
            return 0.0
        return round(1000.0 * len(clashes) / heavy, 2)

    # ------------------------------------------------------------------
    # MolProbity-style van-der-Waals overlap clashes
    # ------------------------------------------------------------------

    @staticmethod
    def _radius(atom, radii):
        el = (atom.element or 'C').upper()
        return radii.get(el, 1.70)

    def vdw_clashes(self, radii=None, overlap_cut:float=None,
                    include_hydrogen:bool=True) -> pd.DataFrame:
        '''vdW-overlap clashes for this chain (MolProbity definition).

        Cutoff: overlap_cut=0.40 A (VDW_OVERLAP_CUT) on overlap =
        (r_i + r_j) - d_ij; same-residue, sequential (i, i+1) and SG-SG
        disulfide pairs excluded.
        '''
        radii = radii or self.VDW_RADII
        overlap_cut = self.VDW_OVERLAP_CUT if overlap_cut is None else overlap_cut
        cols = ['atom_pos', 'atom', 'element', 'other_pos', 'other_atom',
                'other_element', 'distance', 'overlap']
        atoms = [a for r in self.standard_residues()
                 for a in r if (a.element or '').upper() in radii
                 and (include_hydrogen or (a.element or '').upper() != 'H')]
        if len(atoms) < 2:
            return pd.DataFrame(columns=cols)
        search_cut = 2.0 * max(radii.values())     # generous spatial cutoff
        ns = PDB.NeighborSearch(atoms)
        seen, clashes = set(), []
        for atom in atoms:
            for other in ns.search(atom.coord, search_cut, level='A'):
                if atom is other:
                    continue
                key = tuple(sorted((atom.get_serial_number(),
                                    other.get_serial_number())))
                if key in seen:
                    continue
                seen.add(key)
                pos = atom.get_full_id()[3][1]
                other_pos = other.get_full_id()[3][1]
                if pos == other_pos:
                    continue                        # 1-2 / 1-3 bonds
                if abs(pos - other_pos) == 1:
                    continue                        # peptide-bond neighbours
                if atom.get_name() == 'SG' and other.get_name() == 'SG':
                    continue                        # disulfide
                d = float(np.linalg.norm(atom.coord - other.coord))
                overlap = self._radius(atom, radii) \
                    + self._radius(other, radii) - d
                if overlap >= overlap_cut:
                    clashes.append([
                        pos,
                        atom.get_full_id()[4][0],
                        (atom.element or '').upper(),
                        other_pos,
                        other.get_full_id()[4][0],
                        (other.element or '').upper(), d, overlap
                    ])
        return pd.DataFrame(clashes, columns=cols)

    def vdw_clash_score(self, clashes:pd.DataFrame,
                        n_atoms:int=None) -> dict:
        '''clashes per 1000 atoms and overlap stats (MolProbity clashscore).

        Cutoff: clashscore bands <5 very good / 5-10 acceptable / 10-20
        investigate / >20 poor; 0 means no serious overlaps found.
        '''
        heavy = len(self.heavy_atoms()) if n_atoms is None else n_atoms
        if heavy == 0:
            return {'clash_score': 0.0, 'n_clashes': 0, 'max_overlap': 0.0}
        return {
            'clash_score': round(1000.0 * len(clashes) / heavy, 2),
            'n_clashes': len(clashes),
            'max_overlap': round(float(clashes['overlap'].max()), 3)
            if not clashes.empty else 0.0,
        }

    # ------------------------------------------------------------------
    # summary + gates (used by pipelines/qc/chain_clashes.py)
    # ------------------------------------------------------------------

    def run(self, include_hydrogen:bool=False,
            keep_details:bool=True) -> dict:
        '''compute all clash metrics of one chain.

        The steric half is the same metric as Geometry.run() and the
        former backbone_qc_screen(): non-bonded heavy-atom contacts below
        clash_cutoff (2.4 A), severe contacts below
        severe_clash_distance (2.0 A) and the clash score per 1000 heavy
        atoms. The vdW half is the MolProbity clashscore of the former
        structure_qc_screen(): overlaps >= VDW_OVERLAP_CUT (0.40 A) per
        1000 heavy atoms; hydrogens are excluded unless
        include_hydrogen=True (structures with added hydrogens).

        Returns a dict of summary statistics (n_heavy_atoms, num_clashes,
        num_severe_clashes, clash_score, n_vdw_clashes, vdw_clash_score,
        max_overlap); raw tables are nested under 'details' when
        keep_details.
        '''
        heavy = len(self.heavy_atoms())
        steric = self.steric_clashes()
        severe = self.severe_clashes(steric)
        info = {
            'n_heavy_atoms': heavy,
            'num_clashes': len(steric),
            'num_severe_clashes': len(severe),
            'clash_score': self.clash_score(steric, n_atoms=heavy),
        }
        vdw = self.vdw_clashes(include_hydrogen=include_hydrogen)
        score = self.vdw_clash_score(vdw, n_atoms=heavy)
        info.update({
            'n_vdw_clashes': score['n_clashes'],
            'vdw_clash_score': score['clash_score'],
            'max_overlap': score['max_overlap'],
        })
        if keep_details:
            info['details'] = {
                'steric_clashes': steric,
                'vdw_clashes': vdw,
            }
        return info

    def check(self, info:dict=None) -> tuple:
        '''gate on the summary statistics returned by run().

        Cutoffs (fail -> reason):
          num_severe_clashes > 0 -> 'severe_clashes'
            (contacts < severe_clash_distance = 2.0 A);
          clash_score > clash_score_cutoff (30.0/1000 atoms) -> 'clashes';
          vdw_clash_score > vdw_clash_cutoff (20.0/1000 atoms)
            -> 'vdw_clashes' (MolProbity poor band).
        '''
        if info is None:
            info = self.run(keep_details=False)
        if not info:
            return False, 'null_clashes'
        # severe steric clashes: non-bonded heavy-atom contacts below
        # severe_clash_distance (2.0 A)
        if info.get('num_severe_clashes', 0) > 0:
            return False, 'severe_clashes'
        # steric clashes: clash score = clashes per 1000 heavy atoms,
        # see Clashes.clash_score_cutoff
        if info.get('clash_score', 0) > self.clash_score_cutoff:
            return False, 'clashes'
        # MolProbity vdW clashscore above the poor band
        if info.get('vdw_clash_score', 0) > self.vdw_clash_cutoff:
            return False, 'vdw_clashes'
        return True, 'pass'
