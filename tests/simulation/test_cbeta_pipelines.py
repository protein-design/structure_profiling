'''
tests for pipelines/qc/chain_cbeta.py (table chain_cbeta):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the metric
collection itself runs on the local PDB test data (the misplaced-CB
model is synthesized with the PDB-editing helpers in tests/helper.py,
since all shipped PDB test data is clean).
- pytest -s tests/simulation/test_cbeta_pipelines.py
'''
import os
import tempfile
from unittest import TestCase

import numpy as np

from ..helper import *
from ..helper import _atom_name, _is_atom, _resseq
from src.structure_profiling import Pdb, TABLE_COLUMNS

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')


def misplaced_cb_pdb(tmp_dir):
    '''6BSX with the CB of residue 50 (ASN) pulled 1.5 A radially out
    along its own CA->CB direction: cb_dev >> CB_DEVIATION_CUT (0.25 A).'''
    lines = open(PDB_6BSX).readlines()
    coords = {}
    for l in lines:
        if _is_atom(l) and _resseq(l) == 50 and _atom_name(l) in ('CA', 'CB'):
            coords[_atom_name(l)] = pdb_coord(l)
    ca, cb = coords['CA'], coords['CB']
    new_cb = cb + 1.5 * (cb - ca) / np.linalg.norm(cb - ca)
    return write_pdb_lines(
        move_atom(lines, 50, 'CB', new_cb),
        os.path.join(tmp_dir, 'misplaced_cb.pdb'))


class TestCbetaMetrics(TestCase):
    '''cbeta_deviation() summaries via the pipeline (no DB required)'''

    def test_cbeta_metrics_clean(self):
        from pipelines.qc.chain_cbeta import cbeta_metrics
        info = cbeta_metrics(PDB_6BSX, 'A')
        self.assertEqual(info['n_residues'], 170)
        # 10 glycin residues have no CB and cannot deviate
        self.assertEqual(info['n_cb_analysed'], 160)
        self.assertEqual(info['cb_max_dev'], 0.228)
        self.assertEqual(info['cb_outlier_frac'], 0.0)

    def test_cbeta_metrics_misplaced_cb(self):
        from pipelines.qc.chain_cbeta import cbeta_metrics
        with tempfile.TemporaryDirectory() as tmp_dir:
            info = cbeta_metrics(misplaced_cb_pdb(tmp_dir), 'A')
        self.assertEqual(info['n_residues'], 170)
        self.assertEqual(info['n_cb_analysed'], 160)
        self.assertEqual(info['cb_max_dev'], 1.494)
        self.assertEqual(info['cb_outlier_frac'], 0.006)


class TestCbetaPipelineRun(TestCase):
    '''run() wraps cbeta_deviation() into chain_cbeta records'''

    def test_run_yields_table_columns(self):
        from pipelines.qc.chain_cbeta import run
        params = {'table_name': 'chain_cbeta', 'verbose': 0}
        rows = [{
            'pdb_id': '6BSX', 'chain_id': '6BSX_A', 'chain_no': 'A',
            'relative_pdb': '6BSX_0_A_0.pdb', 'pdb_file': PDB_6BSX,
        }]
        records = list(run(rows, params, {}))
        self.assertEqual(len(records), 1)
        rec = dict(zip(TABLE_COLUMNS['chain_cbeta'], records[0]))
        self.assertEqual(rec['pdb_id'], '6BSX')
        self.assertEqual(rec['chain_id'], '6BSX_A')
        self.assertEqual(rec['n_residues'], 170)
        self.assertEqual(rec['n_cb_analysed'], 160)
        self.assertEqual(rec['cb_outlier_frac'], 0.0)

    def test_run_skips_unparsable_files(self):
        from pipelines.qc.chain_cbeta import run
        params = {'table_name': 'chain_cbeta', 'verbose': 0}
        rows = [{'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
                 'relative_pdb': 'nope.pdb',
                 'pdb_file': os.path.join(DATA_DIR, 'nope.pdb')}]
        self.assertEqual(list(run(rows, params, {})), [])


class TestCbetaPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_chain_cbeta(self):
        self.assertEqual(TABLE_COLUMNS['chain_cbeta'],
                         ['pdb_id', 'chain_id', 'n_residues',
                          'n_cb_analysed', 'cb_max_dev', 'cb_outlier_frac'])

    def test_main_pool_registers_chain_cbeta(self):
        from src.app.main import pool
        self.assertEqual(pool['chain_cbeta'], 'qc.chain_cbeta')

    def test_sql_schema_defines_chain_cbeta(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'qc.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE chain_cbeta', sql)
        for col in ('n_residues', 'n_cb_analysed', 'cb_max_dev',
                    'cb_outlier_frac'):
            self.assertIn(col, sql)
