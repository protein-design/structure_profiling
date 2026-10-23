import locale
import os
from pathlib import Path

import pandas as pd
from Bio.PDB.DSSP import DSSP

from ..geometry.pdb import Pdb

DSSP_KEYS = ['dssp_index', 'residue_name', 'dssp_code', 'Relative_ASA', 
    'Phi', 'Psi', 'NH–O_1_relidx', 'NH–O_1_energy', 'O–NH_1_relidx', 
    'O–NH_1_energy', 'NH–O_2_relidx', 'NH–O_2_energy', 'O–NH_2_relidx',
    'O–NH_2_energy', ]

DSSP_CODES ={
    'H': 'alpha helix',
    'B': 'isolated beta bridge',
    'E': 'beta sheet',
    'G': '310 helix',
    'I': 'Pi helix',
    'T': 'hydrogen-bonded turn',
    'S': 'bend', 
    'C': 'coil or unstructured region',
    '-': 'loop',
} 

# mkdssp 4.x (libcifpp based) wants the Chemical Component Dictionary
# (components.cif[.gz]) at runtime: without it every residue becomes an
# "unknown compound" and mkdssp spams a Configuration error on stderr.
# Builds installed from source often have no data dir compiled in, so point
# LIBCIFPP_DATA_DIR at the first location that actually holds the file;
# the mkdssp subprocess spawned by Biopython inherits the environment.
CCD_NAMES = ('components.cif', 'components.cif.gz')
CCD_DIRS = (
    Path.home() / '.local/share/libcifpp',
    Path('/var/cache/libcifpp'),
    Path('/usr/share/libcifpp'),
    Path('/usr/local/share/libcifpp'),
)


class AnalyzeDssp(Pdb):

    def __init__(self, pdb_file:Path, outdir:str=None, verbose:bool=False):
        super().__init__(pdb_file, outdir, verbose)
        self.load_structure()

    @staticmethod
    def ensure_utf8_ctype():
        '''
        native libraries loaded earlier in the same process (e.g. OpenMM)
        may switch the C locale to plain "C"/POSIX. Biopython then decodes
        the mkdssp output as ascii and crashes with UnicodeDecodeError on
        the non-ascii bytes emitted by dssp 4.x. Restore a UTF-8 LC_CTYPE
        when the current one cannot handle them.
        '''
        try:
            current = locale.setlocale(locale.LC_CTYPE, None)
            if current.lower() not in ('c', 'posix', 'ansi_x3.4-1968'):
                return
            for candidate in ('C.UTF-8', 'en_US.UTF-8', 'UTF-8'):
                try:
                    locale.setlocale(locale.LC_CTYPE, candidate)
                    return
                except locale.Error:
                    continue
        except Exception:
            pass

    @staticmethod
    def ensure_ccd_data_dir():
        '''
        set LIBCIFPP_DATA_DIR to the first directory that actually contains
        a CCD file, so the mkdssp subprocess can resolve residue compounds.
        no-op when the variable is already set by the environment.
        '''
        if os.getenv('LIBCIFPP_DATA_DIR'):
            return
        for data_dir in CCD_DIRS:
            if any((Path(data_dir) / name).is_file() for name in CCD_NAMES):
                os.environ['LIBCIFPP_DATA_DIR'] = str(data_dir)
                return

    def dssp_iter(self):
        for model in self.structure:
            try:
                self.ensure_utf8_ctype()
                self.ensure_ccd_data_dir()
                dssp = DSSP(model, self.pdb_file)
                dssp_data = self.parse_dssp(dssp)
                yield model, dssp_data
            except Exception as e:
                print(f"ERROR: DSSP failed for {self.pdb_file}. error={e}")

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