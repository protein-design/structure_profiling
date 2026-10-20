'''
tests for pipelines/qc/chain_chirality.py (table chain_chirality):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the metric
collection itself runs on the local PDB test data (D / planar CA
centres are synthesized with the PDB-editing helpers in
tests/helper.py, since all shipped PDB test data is clean).
- pytest -s tests/simulation/test_chirality_pipelines.py
'''
import os
import tempfile
from unittest import TestCase

import numpy as np

from ..helper import *
from ..helper import _atom_name, _is_atom, _resseq
from src.structure_profiling import Pdb, TABLE_COLUMNS

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')


def _residue_atom(lines, resseq, *names):
    '''coords of the named atoms of one residue.'''
    coords = {}
    for l in lines:
        if _is_atom(l) and _resseq(l) == resseq and _atom_name(l) in names:
            coords[_atom_name(l)] = pdb_coord(l)
    return coords


def d_center_pdb(tmp_dir, resseq=50):
    '''6BSX with the CB of one residue mirrored across its N-CA-C plane:
    the chiral volume flips sign -> a D centre in an L-amino-acid chain
    (CA-CB bond length and all bond angles preserved).'''
    lines = open(PDB_6BSX).readlines()
    co = _residue_atom(lines, resseq, 'N', 'CA', 'C', 'CB')
    N, CA, C, CB = co['N'], co['CA'], co['C'], co['CB']
    normal = np.cross(N - CA, C - CA)
    normal = normal / np.linalg.norm(normal)
    cb_new = CB - 2 * np.dot(CB - CA, normal) * normal
    return write_pdb_lines(
        [set_pdb_coord(l, cb_new)
         if _is_atom(l) and _resseq(l) == resseq and _atom_name(l) == 'CB'
         else l for l in lines],
        os.path.join(tmp_dir, 'd_center.pdb'))


def planar_center_pdb(tmp_dir, resseq=50):
    '''6BSX with the CB of one residue moved into its N-CA-C plane:
    the chiral volume collapses below CHIRAL_VOL_MIN (0.7 A^3) ->
    a 'planar' CA centre.'''
    lines = open(PDB_6BSX).readlines()
    co = _residue_atom(lines, resseq, 'N', 'CA', 'C')
    N, CA, C = co['N'], co['CA'], co['C']
    e1 = (N - CA) / np.linalg.norm(N - CA)
    e2 = (C - CA) / np.linalg.norm(C - CA)
    return write_pdb_lines(
        [set_pdb_coord(l, CA + e1 + e2)
         if _is_atom(l) and _resseq(l) == resseq and _atom_name(l) == 'CB'
         else l for l in lines],
        os.path.join(tmp_dir, 'planar_center.pdb'))


class TestChiralityMetrics(TestCase):
    '''ca_chirality() summaries via the pipeline (no DB required)'''

    def test_chirality_metrics_clean(self):
        from pipelines.qc.chain_chirality import chirality_metrics
        info = chirality_metrics(PDB_6BSX, 'A')
        self.assertEqual(info['n_residues'], 170)
        # 10 glycin residues have no CB and are not analysed ('NA')
        self.assertEqual(info['n_chirality_analysed'], 160)
        self.assertEqual(info['n_chirality_outliers'], 0)
        self.assertEqual(info['chirality_outlier_frac'], 0.0)

    def test_chirality_metrics_d_center(self):
        from pipelines.qc.chain_chirality import chirality_metrics
        with tempfile.TemporaryDirectory() as tmp_dir:
            info = chirality_metrics(d_center_pdb(tmp_dir), 'A')
        self.assertEqual(info['n_chirality_analysed'], 160)
        self.assertEqual(info['n_chirality_outliers'], 1)
        self.assertEqual(info['chirality_outlier_frac'], 0.006)

    def test_chirality_metrics_planar_center(self):
        from pipelines.qc.chain_chirality import chirality_metrics
        with tempfile.TemporaryDirectory() as tmp_dir:
            info = chirality_metrics(planar_center_pdb(tmp_dir), 'A')
        self.assertEqual(info['n_chirality_analysed'], 160)
        self.assertEqual(info['n_chirality_outliers'], 1)
        self.assertEqual(info['chirality_outlier_frac'], 0.006)


class TestChiralityPipelineRun(TestCase):
    '''run() wraps ca_chirality() into chain_chirality records'''

    def test_run_yields_table_columns(self):
        from pipelines.qc.chain_chirality import run
        params = {'table_name': 'chain_chirality', 'verbose': 0}
        rows = [{
            'pdb_id': '6BSX', 'chain_id': '6BSX_A', 'chain_no': 'A',
            'relative_pdb': '6BSX_0_A_0.pdb', 'pdb_file': PDB_6BSX,
        }]
        records = list(run(rows, params, {}))
        self.assertEqual(len(records), 1)
        rec = dict(zip(TABLE_COLUMNS['chain_chirality'], records[0]))
        self.assertEqual(rec['pdb_id'], '6BSX')
        self.assertEqual(rec['chain_id'], '6BSX_A')
        self.assertEqual(rec['n_residues'], 170)
        self.assertEqual(rec['n_chirality_analysed'], 160)
        self.assertEqual(rec['n_chirality_outliers'], 0)

    def test_run_skips_unparsable_files(self):
        from pipelines.qc.chain_chirality import run
        params = {'table_name': 'chain_chirality', 'verbose': 0}
        rows = [{'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
                 'relative_pdb': 'nope.pdb',
                 'pdb_file': os.path.join(DATA_DIR, 'nope.pdb')}]
        self.assertEqual(list(run(rows, params, {})), [])


class TestChiralityPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_chain_chirality(self):
        self.assertEqual(TABLE_COLUMNS['chain_chirality'],
                         ['pdb_id', 'chain_id', 'n_residues',
                          'n_chirality_analysed', 'n_chirality_outliers',
                          'chirality_outlier_frac'])

    def test_main_pool_registers_chain_chirality(self):
        from src.app.main import pool
        self.assertEqual(pool['chain_chirality'], 'qc.chain_chirality')

    def test_sql_schema_defines_chain_chirality(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'qc.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE chain_chirality', sql)
        for col in ('n_residues', 'n_chirality_analysed',
                    'n_chirality_outliers', 'chirality_outlier_frac'):
            self.assertIn(col, sql)

