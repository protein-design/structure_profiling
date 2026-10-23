'''
Ligand / binding-site neighborhood search for one Bio.PDB structure.

Detects the protein residues close in space to a site of interest: a
HETATM ligand (default) or an explicit set of site residues such as an
active site (site_residues). Neighbor residues may or may not be
continuous in the chain sequence - pocket residues are usually
scattered over the sequence - so run() also reports the number of
sequence-continuous fragments they form (n_fragments).

search() returns the per-residue contact table (chain_no, resseq,
resname, min_dist, n_contacts, closest atom pair) exported by
pipelines/dynamics/chain_ligand.py to csv / pickle next to the chain
PDB; run() returns the summary exported to the SQL table chain_ligand.

Cutoffs:
    search_cutoff = 5.0 A   heavy-atom distance flagged as a contact
                            (typical ligand-neighborhood radii 4-6 A)
Waters (hetflag 'W', HOH / DOD / WAT) are skipped as sites unless
include_waters=True.

references:
    Bio.PDB NeighborSearch (KD-tree); PLIP / PDBsum-style
    ligand-neighborhood radii (4-6 A).
'''
import numpy as np
import pandas as pd
from Bio.PDB import NeighborSearch


class Neighbors:
    '''ligand / binding-site neighborhood of one Bio.PDB structure.'''

    # neighborhood cutoff (A); heavy-atom contacts below are neighbors
    search_cutoff:float = 5.0

    # waters are never protein neighbors and are skipped as sites
    WATER_NAMES = ('HOH', 'DOD', 'WAT')

    # columns of the per-residue contact table returned by search()
    # (run() prepends the site label column in front of these)
    NEIGHBOR_COLUMNS = ['chain_no', 'resseq', 'icode', 'resname',
        'min_dist', 'n_contacts', 'atom', 'site_atom']

    def __init__(self, structure, chain_no:str=None):
        '''
        args:
            structure: Bio.PDB.Structure (e.g. Pdb(pdb_file).structure)
            chain_no:  restrict the neighbor candidates to this chain;
                       None keeps every protein chain of the structure
        '''
        self.structure = structure
        self.chain_no = chain_no

    # ------------------------------------------------------------------
    # residue / atom selection
    # ------------------------------------------------------------------

    def protein_residues(self) -> list:
        '''standard amino-acid residues (skip hetero / water), optionally
        restricted to self.chain_no.

        Cutoff: none - filters on residue id[0] == ' ' only (same filter
        as Clashes / Breaks / GeometryDistance).
        '''
        pool = []
        if not self.structure:
            return pool
        for model in self.structure:
            for chain in model:
                if self.chain_no is not None and chain.id != self.chain_no:
                    continue
                pool.extend(r for r in chain if r.id[0] == ' ')
        return pool

    def protein_atoms(self) -> list:
        '''heavy (non-hydrogen) atoms of the standard residues.

        Cutoff: none - filters on element != 'H' only.
        '''
        return [
            a for r in self.protein_residues()
            for a in r if a.element != 'H'
        ]

    def hetero_residues(self, include_waters:bool=False) -> list:
        '''HETATM residues (id[0] != ' '), waters skipped by default.

        Cutoff: none - waters are residues with hetflag 'W' or a water
        resname; they are site candidates only when include_waters.
        Restricted to self.chain_no like the protein candidates: on a
        full complex only the HETATM of that chain become sites.
        '''
        pool = []
        if not self.structure:
            return pool
        for model in self.structure:
            for chain in model:
                if self.chain_no is not None and chain.id != self.chain_no:
                    continue
                for r in chain:
                    if r.id[0] == ' ':
                        continue
                    if not include_waters and self.is_water(r):
                        continue
                    pool.append(r)
        return pool

    def is_water(self, residue) -> bool:
        '''water / solvent residue (hetflag 'W' or water resname).'''
        if residue.id[0] == 'W':
            return True
        return residue.resname.strip().upper() in self.WATER_NAMES

    @staticmethod
    def residue_label(residue) -> str:
        '''compact residue label "chain/resnameresseq[icode]".'''
        icode = residue.id[2] if residue.id[2] else ' '
        label = f"{residue.get_parent().id}/{residue.resname.strip()}{residue.id[1]}"
        return label if icode == ' ' else label + icode.strip()

    def site_residues(self, spec:list) -> list:
        '''resolve an explicit activity site into (label, atoms) pairs.

        args:
            spec: [(chain_no, resseq), ...] - residues defining the site,
                e.g. catalytic residues. Missing residues are skipped
                with a warning.

        Cutoff: none - exact residue (chain_no, resseq) lookup.
        '''
        sites = []
        for item in spec:
            chain_no, resseq = item[0], item[1]
            res = None
            for model in self.structure:
                for chain in model:
                    if chain.id != chain_no:
                        continue
                    for r in chain:
                        if r.id[1] == resseq:
                            res = r
                            break
            if res is None:
                print(f"Warning: site residue {chain_no}/{resseq} not found.")
                continue
            atoms = [a for a in res if a.element != 'H']
            if atoms:
                sites.append((self.residue_label(res), atoms))
        return sites

    def ligand_sites(self, include_waters:bool=False) -> list:
        '''(label, heavy atoms) per HETATM residue - the default sites.

        Cutoff: none - one site per HETATM residue; waters are skipped
        unless include_waters=True.
        '''
        sites = []
        for res in self.hetero_residues(include_waters):
            atoms = [a for a in res if a.element != 'H']
            if atoms:
                sites.append((self.residue_label(res), atoms))
        return sites

    # ------------------------------------------------------------------
    # neighborhood search
    # ------------------------------------------------------------------

    def neighbor_search(self):
        '''KD-tree index of the protein heavy atoms (built once).'''
        if not hasattr(self, '_ns'):
            self._ns = NeighborSearch(self.protein_atoms())
        return self._ns

    def search(self, site_atoms:list, cutoff:float=None) -> pd.DataFrame:
        '''protein residues with any heavy atom within cutoff of the site.

        Cutoff: cutoff=search_cutoff (5.0 A) heavy-atom distance;
        candidates come from a Bio.PDB NeighborSearch (KD-tree) per site
        atom and the exact distances are recomputed on the candidates.
        Hydrogens are ignored on both sides.

        returns DataFrame(chain_no, resseq, icode, resname,
        min_dist, n_contacts, atom, site_atom) sorted by residue; the
        site label column is prepended by run() when sites are combined.
        '''
        cutoff = self.search_cutoff if cutoff is None else cutoff
        site_atoms = [a for a in site_atoms if a.element != 'H']
        if not self.structure or not site_atoms \
                or not self.protein_atoms():
            return pd.DataFrame(columns=self.NEIGHBOR_COLUMNS)

        # KD-tree candidates per site atom, grouped by residue
        ns = self.neighbor_search()
        cand = {}
        for atom in site_atoms:
            for hit in ns.search(atom.coord, cutoff, level='A'):
                res = hit.get_parent()
                cand.setdefault((res.get_parent().id, res.id), set()).add(hit)

        site_names = [a.name for a in site_atoms]
        site_coords = np.asarray([a.coord for a in site_atoms], dtype=float)
        rows = []
        for (chain_no, res_id), res_atoms in cand.items():
            res_atoms = list(res_atoms)
            res_coords = np.asarray([a.coord for a in res_atoms], dtype=float)
            dist = np.linalg.norm(
                res_coords[:, None, :] - site_coords[None, :, :], axis=-1)
            ri, si = np.unravel_index(np.argmin(dist), dist.shape)
            rows.append({
                'chain_no': chain_no,
                'resseq': res_id[1],
                'icode': (res_id[2] or ' ').strip(),
                'resname': res_atoms[0].get_parent().resname.strip(),
                'min_dist': round(float(dist[ri, si]), 3),
                'n_contacts': int((dist.min(axis=1) <= cutoff).sum()),
                'atom': res_atoms[ri].name,
                'site_atom': site_names[si],
            })
        df = pd.DataFrame(rows, columns=self.NEIGHBOR_COLUMNS)
        return df.sort_values(['chain_no', 'resseq', 'icode']) \
            .reset_index(drop=True)

    def fragments(self, neighbors:pd.DataFrame) -> int:
        '''number of sequence-continuous fragments among the neighbors.

        Pocket residues are scattered over the sequence: several
        fragments instead of one stretch. neighbors must be sorted by
        (chain_no, resseq) - a chain change or a numbering gap of > 1
        starts a new fragment.

        Cutoff: none - counts resseq jumps != 1 within a chain.
        '''
        if neighbors is None or neighbors.empty:
            return 0
        n = 0
        prev = None
        for chain_no, resseq in zip(neighbors['chain_no'],
                                    neighbors['resseq']):
            if prev is None or chain_no != prev[0] or resseq != prev[1] + 1:
                n += 1
            prev = (chain_no, resseq)
        return n
    # ------------------------------------------------------------------
    # summary (used by pipelines/dynamics/chain_ligand.py)
    # ------------------------------------------------------------------

    def run(self, cutoff:float=None, site_residues:list=None,
            include_waters:bool=False, keep_details:bool=True) -> dict:
        '''neighborhood summary of one structure.

        The sites are the HETATM ligands (waters skipped unless
        include_waters=True) or the explicit site_residues
        [(chain_no, resseq), ...] of an activity site. Every standard
        residue with a heavy atom within cutoff (search_cutoff = 5.0 A)
        of any site atom is a neighbor; n_neighbor_residues dedupes
        residues touching several sites and n_fragments counts the
        sequence-continuous segments they form.

        Returns a dict of summary statistics (n_sites, ligand_names,
        n_neighbor_residues, n_fragments, min_dist, mean_dist,
        neighbor_residues); the full contact table is nested under
        'neighbors' when keep_details.
        '''
        cutoff = self.search_cutoff if cutoff is None else cutoff
        if site_residues:
            sites = self.site_residues(site_residues)
        else:
            sites = self.ligand_sites(include_waters=include_waters)

        frames = []
        for label, atoms in sites:
            df = self.search(atoms, cutoff)
            df.insert(0, 'site', label)
            frames.append(df)
        neighbors = pd.concat(frames, ignore_index=True) \
            if frames else pd.DataFrame(columns=self.NEIGHBOR_COLUMNS)

        # a site residue is not its own neighbor (activity sites)
        site_keys = {(res.get_parent().id, res.id[1])
            for _, atoms in sites for res in {a.get_parent() for a in atoms}}
        if site_keys and not neighbors.empty:
            keep = [key not in site_keys for key
                    in zip(neighbors['chain_no'], neighbors['resseq'])]
            neighbors = neighbors[keep].reset_index(drop=True)

        # unique neighbor residues over all sites (pocket view)
        uniq = neighbors.drop_duplicates(['chain_no', 'resseq', 'icode']) \
            if not neighbors.empty else neighbors
        info = {
            'n_sites': len(sites),
            'ligand_names': ','.join(dict.fromkeys(
                atoms[0].get_parent().resname.strip()
                for _, atoms in sites)),
            'n_neighbor_residues': len(uniq),
            'n_fragments': self.fragments(uniq),
            'min_dist': round(float(uniq['min_dist'].min()), 3)
                if not uniq.empty else None,
            'mean_dist': round(float(uniq['min_dist'].mean()), 3)
                if not uniq.empty else None,
            'neighbor_residues': ','.join(
                f"{r.chain_no}/{r.resname}{r.resseq}"
                f"{r.icode if r.icode else ''}:{r.min_dist}"
                for r in uniq.itertuples()),
        }
        if keep_details:
            info['neighbors'] = neighbors
        return info

