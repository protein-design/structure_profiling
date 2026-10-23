'''
pLDDT per residue local distance testing
Ported from bio_omics (bioomics.PDB.calculate_plddt): only the parts
used by the predict pipelines are kept; the plotting helpers requiring
matplotlib/seaborn are not ported.
NB: bio_omics provides c_atom_iter / ca_residue_iter on its Pdb class;
structure_profiling's Pdb does not, so they are defined here.
'''
from pathlib import Path
import numpy as np
import pandas as pd

from ..structure.pdb import Pdb


class CalculatePlddt(Pdb):

    def __init__(self, pdb_file, outdir:Path=None, verbose:bool=False):
        super().__init__(pdb_file, outdir, verbose)

    def c_atom_iter(self):
        # alpha C-atom
        for model, chain, residues in self.ca_residue_iter():
            c_atoms = []
            for residue in residues:
                for atom in residue:
                    if atom.id == 'CA':
                        c_atoms.append(atom)
            yield model, chain, c_atoms

    def ca_residue_iter(self):
        '''
        Get all alpha atoms from monomer
        '''
        for model in self.structure:
            for chain in model:
                ca_residues = []
                for residue in chain:
                    # Check if residue has a alpha atom
                    if "CA" in residue: 
                        ca_residues.append(residue)
                yield model, chain, ca_residues

    def complex_ca_bfactor(self) -> dict:
        '''
        scan all models and chains
        '''
        data = {
            'pdb_id': self.structure_id,
            'pdb_file': str(self.pdb_file),
            'bfactor': [],
        }
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            bfactors = []
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                # b-factor is not assay data but determined by refinement process
                bfactor = float(atom.get_bfactor())
                bfactors.append((residue_number, bfactor))
            if bfactors:
                df = pd.DataFrame(bfactors, columns=['residue_number', 'bfactor'])
                stats = self.bfactor_statistics(df)
                item = {
                    'model': model,
                    'chain': chain,
                    'bfactor': df,
                    'stats':stats,
                }
                data['bfactor'].append(item)
        return data

    def bfactor_statistics(self, df:pd.DataFrame, res_start:int=None, res_end:int=None) -> dict:
        """
        Calculate B-factor statistics
            low bfactor: more stable, ordered, or rigid part of structure
            high bfactor: greater aomic displacement and disorder.
        """
        if res_start:
            df=df[(df['residue_number']>=res_start)]
        if res_end:
            df=df[(df['residue_number']<=res_end)]
        num_ca = len(df)
        if num_ca > 1:
            scores = df['bfactor']
            stats = {
                'mean': np.mean(scores),
                'median': np.median(scores),
                'std': np.std(scores),
                'min': np.min(scores),
                'max': np.max(scores),
            }
            return stats
        return {}
        
    def complex_ca_plddt(self) -> dict:
        '''
        scan all models and chains
        '''
        data = {
            'structure_id': self.structure_id,
            'pdb_file': str(self.pdb_file),
            'pLDDT': [],
        }
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            plddt = []
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                bfactor = float(atom.get_bfactor())
                plddt.append((residue_number, bfactor))
            if plddt:
                plddt = pd.DataFrame(plddt, columns=['residue_number', 'pLDDT'])
                stats = self.plddt_statistics(plddt)
                item = {
                    'model': model,
                    'chain': chain,
                    'plddt': plddt,
                    'stats':stats,
                }
                data['pLDDT'].append(item)
        return data
        
    def monomer_ca_plddt(self) -> pd.DataFrame:
        '''
        only retrieve the first one or that is monomer
        '''
        plddt = []
        ca_iter = self.c_atom_iter()
        for model, chain, atoms in ca_iter:
            for atom in atoms:
                residue_number = int(atom.full_id[3][1])
                bfactor = float(atom.get_bfactor())
                plddt.append((residue_number, bfactor))
            break
        if plddt:
            plddt = pd.DataFrame(plddt, columns=['residue_number', 'pLDDT'])
            return plddt

    def plddt_statistics(self, plddt:pd.DataFrame, res_start:int=None, res_end:int=None) -> dict:
        """
        Calculate comprehensive pLDDT statistics
            Very high (pLDDT > 90)
            High (90 > pLDDT > 70)
            Low (70 > pLDDT > 50)
            Very low (pLDDT < 50)
        """
        if res_start:
            plddt=plddt[(plddt['residue_number']>=res_start)]
        if res_end:
            plddt=plddt[(plddt['residue_number']<=res_end)]
        num_ca = len(plddt)
        if num_ca > 1:
            scores = plddt['pLDDT']
            stats = {
                'mean_plddt': np.mean(scores),
                'median_plddt': np.median(scores),
                'std_plddt': np.std(scores),
                'min_plddt': np.min(scores),
                'max_plddt': np.max(scores),
                'confidence_categories': {
                    'very_high': len(scores[scores >= 90]),
                    'high': len(scores[(scores >= 70) & (scores < 90)]),
                    'low': len(scores[(scores >= 50) & (scores < 70)]),
                    'very_low': len(scores[scores < 50])
                },
                'percentage_high_confidence': len(scores[scores >= 70]) / num_ca * 100
            }
            return stats
        return {}

    @staticmethod
    def confidence(plddt_score):
        pool = [(90, 100, 'Very High'), (70, 90, 'High'),
            (50, 70, 'Low'), (0, 50, 'Very Low'),]
        for min_score, max_score, label in pool:
            if min_score <= plddt_score <= max_score:
                return label

    @staticmethod
    def confidence_percentage(plddt, threshold:int):
        names = ('High', 'Low')
        confidence = {i:0 for i in names}
        for score in plddt['pLDDT']:
            if score >= threshold:
                confidence["High"] += 1
            else:
                confidence["Low"] += 1
        n = len(plddt)
        for k in names:
            confidence[k] *= 100 / n
        return confidence