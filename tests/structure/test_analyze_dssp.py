'''
test the moved AnalyzeDssp on example structures
requires the dssp executable (mkdssp)
'''
import os
import shutil
import tempfile
from pathlib import Path
from unittest import TestCase, mock, skipUnless

from src.structure_profiling import AnalyzeDssp
from src.structure_profiling.structure.analyze_dssp import DSSP_CODES

from ..helper import *

pdb_file = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
chain_id = 'A'

HAS_DSSP = shutil.which('mkdssp') is not None or shutil.which('dssp') is not None


@skipUnless(os.path.isfile(pdb_file), f"missing test data: {pdb_file}")
class TestAnalyzeDssp(TestCase):

    @skipUnless(HAS_DSSP, "dssp executable (mkdssp) not installed")
    def test_dssp_iter(self):
        dssp = AnalyzeDssp(pdb_file)
        res = list(dssp.dssp_iter())
        self.assertGreater(len(res), 0)
        model, dssp_data = res[0]
        self.assertIn(chain_id, dssp_data)
        rec = dssp_data[chain_id]
        self.assertGreater(len(rec), 0)
        self.assertIn('dssp_code', rec[0])
        self.assertIn('secondary_structure', rec[0])
        # one-letter dssp alphabet
        for i in rec:
            self.assertIn(i['dssp_code'], DSSP_CODES)

    def test_ensure_ccd_data_dir(self):
        '''LIBCIFPP_DATA_DIR is only set when a CCD file is actually found'''
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch('src.structure_profiling.structure.analyze_dssp.CCD_DIRS', (tmp,)), \
                mock.patch.dict(os.environ):
            os.environ.pop('LIBCIFPP_DATA_DIR', None)
            AnalyzeDssp.ensure_ccd_data_dir()
            self.assertNotIn('LIBCIFPP_DATA_DIR', os.environ)

            Path(tmp, 'components.cif.gz').write_text('fake ccd', encoding='utf-8')
            AnalyzeDssp.ensure_ccd_data_dir()
            self.assertEqual(os.environ['LIBCIFPP_DATA_DIR'], tmp)

            # an explicit user setting always wins
            os.environ['LIBCIFPP_DATA_DIR'] = '/custom'
            AnalyzeDssp.ensure_ccd_data_dir()
            self.assertEqual(os.environ['LIBCIFPP_DATA_DIR'], '/custom')

    def test_parse_dssp_seq(self):
        dssp_seq = 'HHHEEECCCSSTT'
        ft = AnalyzeDssp.parse_dssp_seq(dssp_seq)
        self.assertFalse(ft.empty)
        # 'W' marks the whole chain, skip it when summing segments
        seg = ft[ft['code'] != 'W']
        self.assertEqual(int(seg['length'].sum()), len(dssp_seq))
        codes = set(seg['code'])
        self.assertIn('H', codes)
        self.assertIn('E', codes)

