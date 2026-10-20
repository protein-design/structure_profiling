'''
tests for pipelines/qc/chain_fixed_breaks.py (table chain_fixed_breaks):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the metric
collection itself runs on the local PDB test data (defects are
synthesized with the PDB-editing helpers in tests/helper.py, since
all shipped PDB test data is clean).
- pytest -s tests/simulation/test_fixed_break_pipelines.py
'''
import os
import tempfile
from unittest import TestCase

from ..helper import *
from src.structure_profiling import Pdb, TABLE_COLUMNS

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')


def broken_pdb(tmp_dir):
    '''6BSX with the C-terminal half translated 20 A: one chain break.'''
    lines = open(PDB_6BSX).readlines()
    return write_pdb_lines(
        translate_residues(lines, 87, (20.0, 0.0, 0.0)),
        os.path.join(tmp_dir, 'broken.pdb'))


class TestFixedBreakMetrics(TestCase):
    '''Breaks.run() summaries via the pipeline (no DB required)'''

    def test_break_metrics_clean(self):
        from pipelines.qc.chain_fixed_breaks import break_metrics
        info = break_metrics(PDB_6BSX, 'A')
        self.assertEqual(info['n_residues'], 170)
        self.assertEqual(info['num_missing'], 0)
        self.assertEqual(info['num_breaks'], 0)
        self.assertEqual(info['num_fragments'], 1)
        self.assertEqual(info['num_strained'], 0)
        self.assertEqual(info['ca_max_dist'], 3.88)
        self.assertEqual(info['pep_max_dist'], 1.36)

    def test_break_metrics_broken(self):
        from pipelines.qc.chain_fixed_breaks import break_metrics
        with tempfile.TemporaryDirectory() as tmp_dir:
            info = break_metrics(broken_pdb(tmp_dir), 'A')
        self.assertEqual(info['num_breaks'], 1)
        self.assertEqual(info['num_fragments'], 2)


class TestFixedBreakPipelineRun(TestCase):
    '''run() wraps Breaks.run() into chain_fixed_breaks records'''

    def test_run_yields_table_columns(self):
        from pipelines.qc.chain_fixed_breaks import run
        params = {'table_name': 'chain_fixed_breaks', 'verbose': 0}
        rows = [{
            'pdb_id': '6BSX', 'chain_id': '6BSX_A', 'chain_no': 'A',
            'fixed_pdb': 'chain_fixed/6BSX/6BSX_A/fixed.pdb',
            'pdb_file': PDB_6BSX,
        }]
        records = list(run(rows, params, {}))
        self.assertEqual(len(records), 1)
        rec = dict(zip(TABLE_COLUMNS['chain_fixed_breaks'], records[0]))
        self.assertEqual(rec['pdb_id'], '6BSX')
        self.assertEqual(rec['chain_id'], '6BSX_A')
        self.assertEqual(rec['n_residues'], 170)
        self.assertEqual(rec['num_missing'], 0)
        self.assertEqual(rec['num_breaks'], 0)
        self.assertEqual(rec['num_fragments'], 1)
        self.assertEqual(rec['ca_max_dist'], 3.88)
        self.assertEqual(rec['pep_max_dist'], 1.36)

    def test_run_skips_unparsable_files(self):
        from pipelines.qc.chain_fixed_breaks import run
        params = {'table_name': 'chain_fixed_breaks', 'verbose': 0}
        rows = [{'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
                 'fixed_pdb': 'nope.pdb',
                 'pdb_file': os.path.join(DATA_DIR, 'nope.pdb')}]
        self.assertEqual(list(run(rows, params, {})), [])


class TestFixedBreakPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_chain_fixed_breaks(self):
        self.assertEqual(TABLE_COLUMNS['chain_fixed_breaks'],
                         ['pdb_id', 'chain_id', 'n_residues',
                          'num_missing', 'num_breaks', 'num_fragments',
                          'num_strained', 'ca_max_dist', 'pep_max_dist'])

    def test_main_pool_registers_chain_fixed_breaks(self):
        from src.app.main import pool
        self.assertEqual(pool['chain_fixed_breaks'], 'qc.chain_fixed_breaks')

    def test_sql_schema_defines_chain_fixed_breaks(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'qc.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE chain_fixed_breaks', sql)
        for col in ('n_residues', 'num_missing', 'num_breaks',
                    'num_fragments', 'num_strained', 'ca_max_dist',
                    'pep_max_dist'):
            self.assertIn(col, sql)
