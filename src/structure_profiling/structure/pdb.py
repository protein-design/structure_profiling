'''
process PDB structure: loading and chain access
Decoupled from bio_omics (bioomics.PDB.pdb / parse_structure):
only the parsing used by structure profiling is kept.
'''
from Bio.PDB import PDBParser, MMCIFParser
import gzip
from pathlib import Path

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)


class ParseStructure:
    def __init__(self):
        self.structure = None
        self.chains = []

    def header_values(self, key:str, default=None):
        try:
            return self.structure.header[key]
        except Exception:
            pass
        return default

    def iter_chain(self):
        for model in self.structure:
            for chain in model:
                yield model, chain

    def get_chain(self, chain_no:str=None, model_no:int=0):
        if chain_no is None:
            _, chain = next(self.iter_chain())
            return chain
        for model, chain in self.iter_chain():
            if model.id == model_no and chain.id == chain_no:
                return chain
        return None


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