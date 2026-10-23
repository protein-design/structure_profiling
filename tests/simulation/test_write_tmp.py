'''
tests for structure_profiling/common/write_tmp.py (WriteTmp):
the af_msa pipeline runs colabfold in a ramdisk temp dir; the results
must only be moved to the local outdir when the run succeeded
(move_temp defaults to True). A failed run (move_temp=False) must not
copy anything to the destination - the temp results are discarded.
- pytest -s tests/simulation/test_write_tmp.py
'''
import os
import shutil
import tempfile
from unittest import TestCase

from src.structure_profiling import WriteTmp


class TestWriteTmp(TestCase):

    def setUp(self):
        self.tmp_root = tempfile.mkdtemp(prefix='write_tmp_test_')
        self.destination = os.path.join(self.tmp_root, 'destination')

    def tearDown(self):
        shutil.rmtree(self.tmp_root, ignore_errors=True)

    def test_default_moves_results(self):
        # default behaviour: results are moved to the destination on exit
        with WriteTmp(self.tmp_root, self.destination, 'io') as writer:
            temp_file = writer.temp_path / 'result.a3m'
            temp_file.write_text('>query\nMKVL\n')
            self.assertTrue(temp_file.is_file())
        self.assertTrue(
            os.path.isfile(os.path.join(self.destination, 'result.a3m')))

    def test_skip_move_discards_failed_results(self):
        # af_msa failed run: nothing is copied to the destination,
        # the temporary results are discarded with the temp dir
        with WriteTmp(self.tmp_root, self.destination, 'io') as writer:
            temp_file = writer.temp_path / 'result.a3m'
            temp_file.write_text('>partial\nMK\n')
            writer.move_temp = False
        self.assertFalse(os.path.exists(
            os.path.join(self.destination, 'result.a3m')))

    def test_destination_dir_created(self):
        # the destination dir is created eagerly in the constructor
        destination = os.path.join(self.tmp_root, 'nested', 'out')
        WriteTmp(self.tmp_root, destination, 'io')
        self.assertTrue(os.path.isdir(destination))
