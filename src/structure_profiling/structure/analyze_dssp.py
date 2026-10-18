import pandas as pd
from Bio.PDB.DSSP import DSSP
from .pdb import Pdb

DSSP_KEYS = ['dssp_index', 'residue_name', 'dssp_code', 'Relative_ASA', 
    'Phi', 'Psi', 'NH–O_1_relidx', 'NH–O_1_energy', 'O–NH_1_relidx', 
    'O–NH_1_energy', 'NH–O_2_relidx', 'NH–O_2_energy', 'O–NH_2_relidx',
    'O–NH_2_energy', ]

DSSP_CODES ={'H': 'alpha helix', 'B': 'isolated beta bridge',
    'E': 'beta sheet', 'G': '310 helix', 'I': 'Pi helix',
    'T': 'hydrogen-bonded turn', 'S': 'bend', 
    'C': 'coil or unstructured region',
    '-': 'loop',} 

class AnalyzeDssp(Pdb):

    def __init__(self, pdb_file:str, outdir:str=None, verbose:bool=False):
        super().__init__(pdb_file, outdir, verbose)
        self.load_structure()
    
    def dssp_iter(self):
        for model in self.structure:
            try:
                dssp = DSSP(model, self.pdb_file)
                dssp_data = self.parse_dssp(dssp)
                yield model, dssp_data
            except Exception as e:
                pass

    def parse_dssp(self, dssp):
        res = {}
        for key in dssp.keys():
            chain_id = key[0]
            if chain_id not in res:
                res[chain_id] = []
            values = dssp[key]
            item = {k:v for k,v in zip(DSSP_KEYS, values)}
            item['secondary_structure'] = DSSP_CODES.get(item['dssp_code'])
            item['chain_id'] = chain_id
            item['residue_id'] = key[1][1]
            res[chain_id].append(item)
        return res

    @staticmethod
    def parse_dssp_seq(dssp_seq):
        idx = 0
        last_code = dssp_seq[idx]
        res = {
            'W': [[0, len(dssp_seq)],],
            last_code: [[idx,],],
        }
        for i in range(idx, len(dssp_seq)):
            if dssp_seq[i] == last_code:
                res[last_code][-1].append(i)
            else:
                last_code = dssp_seq[i]
                if last_code not in res:
                    res[last_code] = []
                res[last_code].append([i,])
        # format
        ft = []
        for dssp_code, pos_pool in res.items():
            for pos in pos_pool:
                item = {
                    'code': dssp_code,
                    'structure': DSSP_CODES.get(dssp_code),
                    'start': pos[0],
                    'end': pos[-1] + 1,
                    'length': pos[-1] - pos[0] + 1,
                }
                if dssp_code == 'w':
                    item['structure'] = 'whole chain'
                ft.append(item)
        ft = pd.DataFrame(ft)
        return ft