'''
Parse PDB structure
only the parsing used by structure profiling is kept.
'''
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
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
