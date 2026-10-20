'''
process PDB structure: loading and chain access
Decoupled from bio_omics (bioomics.PDB.pdb / parse_structure):
only the parsing used by structure profiling is kept.
'''
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
import gzip
from pathlib import Path

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

from .atom_select import ChainAtomSelect, LigandAtomSelect, AllAtomSelect

from .parse_structure import ParseStructure

class Pdb(ParseStructure):

    def __init__(self, pdb_file:Path, outdir:Path=None, verbose:bool=False):
        super().__init__()
        # pdb_file, accept both str and Path
        self.pdb_file = Path(pdb_file)
        self.outdir = Path(outdir) if outdir else self.pdb_file.parent
        self.verbose = verbose
        # load structure
        self.load_structure()

    def load_pdb(self):
        self.parser = PDBParser(PERMISSIVE=1)
        self.structure = self.parser.get_structure(
            self.pdb_file.name, str(self.pdb_file)
        )

    def load_pdb_gz(self):
        self.parser = PDBParser(PERMISSIVE=1)
        with gzip.open(self.pdb_file, "rt") as f:
            self.structure = self.parser.get_structure(self.pdb_file.name, f)

    def load_cif(self):
        self.parser = MMCIFParser(QUIET=True, auth_chains=True)
        self.structure = self.parser.get_structure(
            self.pdb_file.name, str(self.pdb_file)
        )

    def load_cif_gz(self):
        self.parser = MMCIFParser(QUIET=True, auth_chains=False)
        with gzip.open(self.pdb_file, "rt") as f:
             self.structure = self.parser.get_structure(self.pdb_file.name, f)

    def load_structure(self):
        if not self.pdb_file.is_file():
            return None
        try:
            if self.pdb_file.name.endswith('.cif.gz'):
                self.load_cif_gz()
            elif self.pdb_file.name.endswith('.cif'):
                self.load_cif()
            elif self.pdb_file.name.endswith('.pdb.gz') or self.pdb_file.name.endswith('.ent.gz'):
                self.load_pdb_gz()
            else:
                self.load_pdb()

            if self.structure:
                self.structure_id = self.structure.header.get('idcode', '').upper()
                if self.verbose:
                    print(f"Successfully retrieve structure of {self.structure_id}")
        except Exception as e:
            print(f"Error: {self.pdb_file} is invalid. error={e}")
            self.structure = []
            self.structure_id = None

    def header_lines(self):
        header_lines = []
        if self.pdb_file.suffix == '.gz':
            with gzip.open(self.pdb_file, 'rt') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)
        else:
            with open(self.pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)
        return header_lines

    def to_chain_pdb(self, model_no:int, chain_no_pool:list, outfile):
        """
        Save chain while preserving original header information
        Decoupled from bio_omics (bioomics.PDB.pdb) for the combo2..5 pipelines.
        """
        try:
            # Write output with header
            io = PDBIO(use_model_flag=model_no)
            io.set_structure(self.structure)
            # selector
            _ids = ''.join(chain_no_pool)
            atom_selector = ChainAtomSelect
            if '-' in _ids:
                atom_selector = LigandAtomSelect
            elif '+' in _ids:
                atom_selector = AllAtomSelect

            # Copy header from original file
            header_lines = self.header_lines()
            with open(outfile, "w") as f:
                f.writelines(header_lines)
                io.save(f, atom_selector(chain_no_pool))
            if self.verbose:
                print(f"Saved chain {chain_no_pool} to {outfile}")
            return True
        except Exception as e:
            print(f'Error: can not extract chain {chain_no_pool} from {self.structure_id}, error={e}')
        return False