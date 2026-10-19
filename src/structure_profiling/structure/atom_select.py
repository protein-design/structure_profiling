'''
select atoms of given chains while exporting a PDB file.
Decoupled from bio_omics (bioomics.PDB.atom_select): only the selectors
used by Pdb.to_chain_pdb are kept.
'''
from Bio.PDB import Select


class ChainAtomSelect(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return True
        return False

    def accept_residue(self, residue):
        if residue.get_id()[0] == ' ':
            return True
        return False


class LigandAtomSelect(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return True
        return False

    def accept_residue(self, residue):
        if residue.get_id()[0] != ' ':
            return True
        return False


class AllAtomSelect(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return True
        return False

    def accept_residue(self, residue):
        if residue.get_id():
            return True
        return False