'''
evaluate predicted structures against sample structures
Decoupled from bio_omics (bioomics.simulation.eval_pred):
- batch_rmsd uses the Biopython Superimposer instead of pyrosetta
- TMalign keeps the MMalign binary interface
prodigy is not kept (stays in bio_omics for the design pipelines).
'''
import subprocess
from pathlib import Path

from Bio.PDB import PDBParser, Superimposer


class EvalPred:

    def __init__(self, bin_dir:Path=None):
        self.bin_dir = bin_dir

    @staticmethod
    def _ca_atom_map(pdb_file) -> dict:
        '''
        map (residue_number, insertion_code) -> CA atom
        '''
        parser = PDBParser(QUIET=True, PERMISSIVE=1)
        structure = parser.get_structure(
            Path(pdb_file).name, str(pdb_file)
        )
        model = next(structure.get_models())
        atoms = {}
        for residue in model.get_residues():
            if 'CA' in residue:
                hetflag, resseq, icode = residue.id
                atoms[(resseq, icode)] = residue['CA']
        return atoms

    def rmsd(self, native_pdb:str, model_pdb:str):
        '''
        CA RMSD after optimal superposition
        '''
        try:
            native = self._ca_atom_map(native_pdb)
            model = self._ca_atom_map(model_pdb)
            keys = sorted(set(native) & set(model))
            if len(keys) < 2:
                return None
            sup = Superimposer()
            sup.set_atoms([native[k] for k in keys],
                [model[k] for k in keys])
            return sup.rms
        except Exception as e:
            print(f"ERROR: rmsd of {native_pdb} vs {model_pdb}, error={e}")
        return None

    def batch_rmsd(self, pool:list):
        '''
        yield each row of pool with row['rmsd'] added,
        pairing residues by residue number
        '''
        for row in pool:
            rmsd = self.rmsd(row['pdb1'], row['pdb2'])
            if rmsd is not None:
                row['rmsd'] = rmsd
                yield row
    
    def TMalign(self, pdb1:str, pdb2:str):
        try:
            exe = str(self.bin_dir / 'MMalign')
            cmd = [exe, str(pdb1), str(pdb2), '-outfmt', '2',]
            result = subprocess.run(cmd, capture_output=True,
                text=True, check=True)
            lines = result.stdout.split('\n')
            l0 = lines[0].split('\t')
            l1 = lines[1].split('\t')
            res = {k:v for k,v in zip(l0, l1)}
            return res
        except Exception as e:
            print(e)
            pass
        return {}