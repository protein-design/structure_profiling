'''
Fix PDB raw data for MD simulation
PDBFixer: https://github.com/openmm/pdbfixer
'''
from datetime import datetime
import gzip
from pathlib import Path
from pdbfixer import PDBFixer
from openmm import app

from .timeout import timeout
from ..common import ProcessJson

class OpenmmPdbfixer:
    def __init__(self, pdb_file:Path, args:dict):
        self.pdb_file = pdb_file
        self.verbose = args.get('verbose', False)
        self.outdir = args.get('outdir')
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.info = {
            'status': 'fail',
            'pdb_file': self.pdb_file,
        }
        self.fixer = None
        self.init_fixer()

    def init_fixer(self):
        try:
            # gzipped PDB: Path.stem of 'x.pdb.gz' is 'x.pdb', so test the
            # suffix instead
            if str(self.pdb_file).endswith('.gz'):
                with gzip.open(self.pdb_file, 'rt') as f:
                    self.fixer = PDBFixer(pdbfile=f)
            else:
                with open(self.pdb_file, 'r') as f:
                    self.fixer = PDBFixer(pdbfile=f)
        except Exception as e:
            print(f"Error: wrong pdb file: {self.pdb_file}")
            self.info['error_init_fixer'] = str(e)

    @timeout(10)
    def fast_fixer(self):
        '''
        timeout is 10 minutes for large-scale fixing
        '''
        try:
            if self.verbose:
                print(f"\nTry to fix {self.pdb_file}:")
            start_time = datetime.now()
            
            # 1. fix PDB
            self.fix_pdb()
            # 2. save data
            self.save_fixer()
            self.info['status'] = 'succeed'

            delta = datetime.now() - start_time
            self.info['duration_min'] = round(delta.total_seconds() / 60)
            if self.verbose:
                print(f"Fixing {self.pdb_file} is done | {datetime.now()}")
        except TimeoutError as e:
            self.info['timeout_fix_pdb'] = str(e)
        except Exception as e:
            self.info['fail_fix_pdb'] = str(e)
        self.save_info()
        return self.info

    def fix_pdb(self):
        if self.fixer is None:
            return None
        self.info.update({
            'num_chains': len(list(self.fixer.topology.chains())),
            "num_residues": len(list(self.fixer.topology.residues())),
            "num_atoms": len(list(self.fixer.topology.atoms())),
        })
        # Identify missing residues
        self.missing_residues()
        # Find non-standard residues
        self.non_standard_residues()
        # remove heterogens, Remove water and ions (keep=False)
        self.remove_heterogens()
        # Find missing atoms
        self.missing_atoms()
        # Add missing hydrogens at pH 7.0
        self.add_hydrogens()

    def missing_residues(self):
        try:
            if self.verbose:
                print(f"\tFix missing residues {datetime.now()}")
            self.fixer.findMissingResidues()
            self.info['missing_residues'] = len(self.fixer.missingResidues) \
                if self.fixer.missingResidues else 0
        except Exception as e:
            self.info['error_missing_residues'] = str(e)

    def remove_heterogens(self):
        try:
            if self.verbose:
                print(f"\tRemove heterogens {datetime.now()}")
            self.fixer.removeHeterogens(keepWater=False)
        except Exception as e:
            self.info['error_remove_hetrogens'] = str(e)

    def non_standard_residues(self):
        try:
            if self.verbose:
                print(f"\tFix non standard residues {datetime.now()}")
            self.fixer.findNonstandardResidues()
            self.info['none_starndard_residues'] = len(self.fixer.nonstandardResidues) \
                if self.fixer.nonstandardResidues else 0
            # Replace them with standard ones
            self.fixer.replaceNonstandardResidues()
        except Exception as e:
            self.info['error_non_standard_residues'] = str(e)

    def missing_atoms(self):
        try:
            if self.verbose:
                print(f"\tFix missing atoms {datetime.now()}")
            self.fixer.findMissingAtoms()
            self.info['missing_atoms'] = len(list(self.fixer.missingAtoms))
            # Add missing atoms
            self.fixer.addMissingAtoms()     
        except Exception as e:
            self.info['error_missing_atoms'] = str(e)

    def add_hydrogens(self):
        try:
            if self.verbose:
                print(f"\tAdd hydrogens {datetime.now()}")
            self.fixer.addMissingHydrogens(7.0)       
            self.info['pH'] = 7.0
        except Exception as e:
            self.info['error_add_hydrogens'] = str(e)

    def save_fixer(self):
        outfile = self.outdir / 'fixed.pdb'
        # Save the fixed structure
        try:
            if self.fixer:
                with open(outfile, 'w') as f:
                    app.PDBFile.writeFile(self.fixer.topology, self.fixer.positions, f)
                self.info['fixed_pdb'] = outfile
        except Exception as e:
            print(f"Warning: cannot save fixed pdb to {outfile}")
            self.info['error_save_fixer'] = str(e)
    
    def save_info(self):
        outfile = self.outdir / 'fixed_info.json'
        ProcessJson.save_json(self.info, outfile)
        if self.verbose:
            print(self.info)
