'''
tests for pipelines/qc/calc_clashes.py (table pdb_clashes):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the metric
collection itself runs on the local PDB test data.
- pytest -s tests/simulation/test_clash_pipelines.py
'''
import os
from unittest import TestCase, skipUnless

from ..helper import *
from src.structure_profiling import TABLE_COLUMNS

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestClashMetrics(TestCase):
    '''Clashes.run() summaries via the pipeline (no DB required)'''

    def test_clash_metrics_strained(self):
        from pipelines.qc.chain_clashes import clash_metrics
        info = clash_metrics(PDB_2OX8, 'A')
        self.assertEqual(info['n_heavy_atoms'], 1078)
        self.assertEqual(info['num_clashes'], 2)
        self.assertEqual(info['num_severe_clashes'], 0)
        self.assertEqual(info['clash_score'], 1.86)
        self.assertEqual(info['n_vdw_clashes'], 25)
        self.assertEqual(info['vdw_clash_score'], 23.19)

    def test_clash_metrics_clean(self):
        from pipelines.qc.chain_clashes import clash_metrics
        info = clash_metrics(PDB_6BSX, 'A')
        self.assertEqual(info['num_clashes'], 0)
        self.assertEqual(info['num_severe_clashes'], 0)
        self.assertEqual(info['vdw_clash_score'], 3.79)


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestClashPipelineRun(TestCase):
    '''run() wraps Clashes.run() into pdb_clashes records'''

    def test_run_yields_table_columns(self):
        from pipelines.qc.chain_clashes import run
        params = {'table_name': 'pdb_clashes', 'host': 'test-host'}
        rows = [{
            'pdb_id': '2OX8', 'chain_id': '2OX8_A', 'chain_no': 'A',
            'relative_pdb': '2OX8_A.pdb', 'pdb_file': PDB_2OX8,
        }]
        records = list(run(rows, params, {}))
        self.assertEqual(len(records), 1)
        rec = dict(zip(TABLE_COLUMNS['pdb_clashes'], records[0]))
        self.assertEqual(rec['pdb_id'], '2OX8')
        self.assertEqual(rec['chain_id'], '2OX8_A')
        self.assertEqual(rec['host'], 'test-host')
        self.assertEqual(rec['n_heavy_atoms'], 1078)
        self.assertEqual(rec['num_clashes'], 2)
        self.assertEqual(rec['vdw_clash_score'], 23.19)

    def test_run_skips_unparsable_files(self):
        from pipelines.qc.chain_clashes import run
        params = {'table_name': 'pdb_clashes', 'host': 'test-host'}
        rows = [{'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
                 'relative_pdb': 'nope.pdb',
                 'pdb_file': os.path.join(DATA_DIR, 'nope.pdb')}]
        self.assertEqual(list(run(rows, params, {})), [])


class TestClashPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_pdb_clashes(self):
        self.assertEqual(TABLE_COLUMNS['pdb_clashes'],
                         ['pdb_id', 'chain_id', 'n_heavy_atoms',
                          'num_clashes', 'num_severe_clashes',
                          'clash_score', 'n_vdw_clashes',
                          'vdw_clash_score', 'max_overlap', 'host'])

    def test_main_pool_registers_calc_clashes(self):
        from src.app.main import pool
        self.assertEqual(pool['calc_clashes'], 'qc.calc_clashes')

    def test_sql_schema_defines_pdb_clashes(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'simulate.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE pdb_clashes', sql)
        for col in ('n_heavy_atoms', 'num_severe_clashes', 'clash_score',
                    'n_vdw_clashes', 'vdw_clash_score', 'max_overlap'):
            self.assertIn(col, sql)
