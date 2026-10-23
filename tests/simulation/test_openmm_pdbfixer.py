'''
integration tests for OpenmmPdbfixer on a small AlphaFold model
(tests/data/lysozyme/AF-P61626-F1-model_v4.pdb)
- pytest -s tests/simulation/test_openmm_pdbfixer.py
'''
import gzip
import os
import shutil
import tempfile
from pathlib import Path

from ..helper import *
from src.structure_profiling.structure.openmm_pdbfixer import OpenmmPdbfixer

PDB_LYSOZYME = os.path.join(
    DATA_DIR, 'lysozyme', 'AF-P61626-F1-model_v4.pdb')


@ddt
class TestOpenmmPdbfixer(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir = Path(tempfile.mkdtemp(prefix='bioomics_fixer_'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir, ignore_errors=True)

    def test_init(self):
        obj = OpenmmPdbfixer(Path(PDB_LYSOZYME), {'outdir': self.outdir})
        self.assertIsNotNone(obj.fixer)
        self.assertEqual(obj.info['status'], 'fail')
        self.assertEqual(obj.info['pdb_file'], Path(PDB_LYSOZYME))

    def test_init_gzipped_pdb(self):
        # regression: the .gz check used Path.stem, which strips '.gz'
        # and made the gzip branch dead code
        gz_file = self.outdir / 'lysozyme.pdb.gz'
        with open(PDB_LYSOZYME, 'rb') as fin, gzip.open(gz_file, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        obj = OpenmmPdbfixer(gz_file, {'outdir': self.outdir})
        self.assertIsNotNone(obj.fixer)

    def test_init_wrong_file(self):
        obj = OpenmmPdbfixer(Path('/not/a/pdb.pdb'), {'outdir': self.outdir})
        self.assertIsNone(obj.fixer)
        self.assertIn('error_init_fixer', obj.info)

    def test_fix_pdb(self):
        obj = OpenmmPdbfixer(Path(PDB_LYSOZYME), {'outdir': self.outdir})
        obj.fix_pdb()
        self.assertGreater(obj.info['num_residues'], 100)
        self.assertGreater(obj.info['num_atoms'], obj.info['num_residues'])
        # AlphaFold models are complete: nothing missing at pH 7.0
        self.assertEqual(obj.info['missing_residues'], 0)
        self.assertEqual(obj.info['pH'], 7.0)

    def test_fast_fixer(self):
        obj = OpenmmPdbfixer(Path(PDB_LYSOZYME), {'outdir': self.outdir})
        info = obj.fast_fixer()
        self.assertEqual(info['status'], 'succeed')
        self.assertTrue((self.outdir / 'fixed.pdb').exists())
        self.assertTrue((self.outdir / 'fixed_info.json').exists())
