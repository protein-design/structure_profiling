'''
tests for the equal_seq_rmsd pipeline moved from bio_omics:
- calculate_rmsd runs without a database (local PDB files only,
  EvalPred uses the Biopython Superimposer)
- TABLE_COLUMNS and the app pool cover the equal_seq_rmsd table
Importing the pipeline module requires a local .env (output_dir etc.)
used by src.helper at import time.
'''
import shutil
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest import skipUnless

project_dir = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_dir))

HAS_ENV = (project_dir / '.env').is_file()
PDB_FILE = project_dir / 'tests' / 'data' / 'PDB' / '2OX8_A.pdb'


@skipUnless(PDB_FILE.is_file(), 'requires tests/data/PDB/2OX8_A.pdb')
class TestCalculateRmsd(unittest.TestCase):
    '''calculate_rmsd pairs chain pdbs and calls EvalPred.batch_rmsd'''

    def setUp(self):
        from ..helper import translate_residues
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp_dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        # chain A: identical copy, chain B: second half of the chain
        # translated (residues are numbered 606..734), so the pair is
        # deformed and superposition cannot cancel the move
        shutil.copy(PDB_FILE, self.tmp_dir / 'a.pdb')
        lines = PDB_FILE.read_text().splitlines(keepends=True)
        moved = self.tmp_dir / 'b.pdb'
        moved.write_text(''.join(translate_residues(lines, 670, (20, 0, 0))))

    def run_calculate_rmsd(self, chain_pdb, ids_pool):
        from src.pipelines.dynamics import equal_seq_rmsd
        params = {'output_pdb_dir': self.tmp_dir}
        meta = defaultdict(int)
        recs = list(equal_seq_rmsd.calculate_rmsd(
            ids_pool, chain_pdb, params, meta))
        return recs, meta

    def test_identical_pair(self):
        recs, meta = self.run_calculate_rmsd(
            {'A': 'a.pdb', 'B': 'a.pdb'}, [('A', 'B')])
        self.assertEqual(len(recs), 1)
        id1, id2, rmsd, rel1, rel2 = recs[0]
        self.assertEqual((id1, id2), ('A', 'B'))
        self.assertAlmostEqual(rmsd, 0.0, places=6)
        self.assertEqual((rel1, rel2), ('a.pdb', 'a.pdb'))
        self.assertEqual(meta['succeed'], 1)
        self.assertEqual(meta['fail'], 0)

    def test_moved_pair(self):
        recs, meta = self.run_calculate_rmsd(
            {'A': 'a.pdb', 'B': 'b.pdb'}, [('A', 'B')])
        self.assertEqual(len(recs), 1)
        self.assertGreater(recs[0][2], 0.1)
        self.assertEqual(meta['succeed'], 1)

    def test_missing_pdbs_are_skipped(self):
        recs, meta = self.run_calculate_rmsd(
            {'A': 'a.pdb'}, [('A', 'B'), ('C', 'D')])
        self.assertEqual(recs, [])
        self.assertEqual(meta['succeed'], 0)
        self.assertEqual(meta['fail'], 0)


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestEqualSeqRmsdRegistry(unittest.TestCase):

    def test_pipeline_import(self):
        from src.pipelines.dynamics import equal_seq_rmsd
        for name in ('validate', 'pull_chain_pdb', 'pull_equal_rmsd',
                'pull_equal_ids', 'calculate_rmsd'):
            self.assertTrue(callable(getattr(equal_seq_rmsd, name)))

    def test_table_columns_cover_equal_seq_rmsd(self):
        from src.structure_profiling.database.build_complex import TABLE_COLUMNS
        self.assertEqual(TABLE_COLUMNS['equal_seq_rmsd'],
            ['chain_id', 'equal_chain_id', 'rmsd',
             'chain_pdb', 'equal_chain_pdb'])

    def test_app_pool_entry(self):
        from src.app import main
        self.assertEqual(main.pool['equal_seq_rmsd'],
            'dynamics.equal_seq_rmsd')