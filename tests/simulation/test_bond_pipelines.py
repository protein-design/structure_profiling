'''
tests for pipelines/qc/chain_bonds.py (table chain_bonds):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the metric
collection itself runs on the local PDB test data (the strained
covalent-geometry model is synthesized with the PDB-editing helpers
in tests/helper.py, since all shipped PDB test data is clean).
- pytest -s tests/simulation/test_bond_pipelines.py
'''
import os
import tempfile
from unittest import TestCase

import numpy as np

from ..helper import *
from src.structure_profiling import Pdb, Geometry, TABLE_COLUMNS

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')


def strained_pdb(tmp_dir):
    '''6BSX with peptide 49-50 stretched by 0.5 A: the C-N bond leaves
    the Engh & Huber window (|Z| >> 5) without breaking the chain.'''
    lines = open(PDB_6BSX).readlines()
    return write_pdb_lines(
        stretch_peptide(lines, 50, extra=0.5),
        os.path.join(tmp_dir, 'strained.pdb'))


class TestBondMetrics(TestCase):
    '''bond-length summaries via the pipeline (no DB required)'''

    def test_bond_metrics_clean(self):
        from pipelines.qc.chain_bonds import bond_metrics
        info = bond_metrics(PDB_6BSX, 'A')
        self.assertEqual(info['n_bonds'], 839)
        self.assertAlmostEqual(info['bond_rmsz'], 0.6415, places=3)
        self.assertEqual(info['n_bond_outliers'], 0)
        self.assertEqual(info['bond_z_outlier_frac'], 0.0)
        self.assertEqual(info['bond_max_absz'], 2.369)

    def test_bond_metrics_rmsz_matches_chain_rmsz(self):
        '''bond_rmsz must equal the kind='bond' half of chain_rmsz()'''
        from pipelines.qc.chain_bonds import bond_metrics
        info = bond_metrics(PDB_6BSX, 'A')
        rmsz = Geometry(PDB_6BSX).chain_rmsz(Pdb(PDB_6BSX).get_chain('A'))
        self.assertAlmostEqual(info['bond_rmsz'], rmsz['bond_rmsz'])
        self.assertEqual(info['n_bond_outliers'], rmsz['n_bond_outliers'])

    def test_bond_metrics_strained(self):
        from pipelines.qc.chain_bonds import bond_metrics
        with tempfile.TemporaryDirectory() as tmp_dir:
            info = bond_metrics(strained_pdb(tmp_dir), 'A')
        self.assertEqual(info['n_bonds'], 839)
        self.assertAlmostEqual(info['bond_rmsz'], 1.464, places=3)
        self.assertEqual(info['n_bond_outliers'], 2)
        self.assertEqual(info['bond_z_outlier_frac'], 0.002)
        self.assertEqual(info['bond_max_absz'], 36.697)

    def test_bond_metrics_clean_2ox8(self):
        from pipelines.qc.chain_bonds import bond_metrics
        info = bond_metrics(PDB_2OX8, 'A')
        self.assertEqual(info['n_bonds'], 634)
        self.assertAlmostEqual(info['bond_rmsz'], 0.383, places=3)
        self.assertEqual(info['n_bond_outliers'], 0)


class TestBondPipelineRun(TestCase):
    '''run() wraps bond_length_z() into chain_bonds records'''

    def test_run_yields_table_columns(self):
        from pipelines.qc.chain_bonds import run
        params = {'table_name': 'chain_bonds', 'verbose': 0}
        rows = [{
            'pdb_id': '6BSX', 'chain_id': '6BSX_A', 'chain_no': 'A',
            'relative_pdb': '6BSX_0_A_0.pdb', 'pdb_file': PDB_6BSX,
        }]
        records = list(run(rows, params, {}))
        self.assertEqual(len(records), 1)
        rec = dict(zip(TABLE_COLUMNS['chain_bonds'], records[0]))
        self.assertEqual(rec['pdb_id'], '6BSX')
        self.assertEqual(rec['chain_id'], '6BSX_A')
        self.assertEqual(rec['n_bonds'], 839)
        self.assertAlmostEqual(rec['bond_rmsz'], 0.6415, places=3)
        self.assertEqual(rec['n_bond_outliers'], 0)

    def test_run_skips_unparsable_files(self):
        from pipelines.qc.chain_bonds import run
        params = {'table_name': 'chain_bonds', 'verbose': 0}
        rows = [{'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
                 'relative_pdb': 'nope.pdb',
                 'pdb_file': os.path.join(DATA_DIR, 'nope.pdb')}]
        self.assertEqual(list(run(rows, params, {})), [])


class TestBondPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_chain_bonds(self):
        self.assertEqual(TABLE_COLUMNS['chain_bonds'],
                         ['pdb_id', 'chain_id', 'n_bonds', 'bond_rmsz',
                          'n_bond_outliers', 'bond_z_outlier_frac',
                          'bond_max_absz'])

    def test_main_pool_registers_chain_bonds(self):
        from src.app.main import pool
        self.assertEqual(pool['chain_bonds'], 'qc.chain_bonds')

    def test_sql_schema_defines_chain_bonds(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'qc.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE chain_bonds', sql)
        for col in ('n_bonds', 'bond_rmsz', 'n_bond_outliers',
                    'bond_z_outlier_frac', 'bond_max_absz'):
            self.assertIn(col, sql)
