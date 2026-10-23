'''
tests for pipelines/dynamics/chain_ligand.py (table chain_ligand):
module import must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time; the neighborhood
collection itself runs on a small synthetic chain PDB with a HETATM
ligand (tests.helper.synthesize_site_pdb, exact coordinates).
- pytest -s tests/simulation/test_ligand_pipelines.py
'''
import os
import pickle
import tempfile
from collections import defaultdict
from pathlib import Path
from unittest import TestCase, skipUnless

import pandas as pd

from ..helper import DATA_DIR, profiling_dir, synthesize_site_pdb
from src.structure_profiling import TABLE_COLUMNS

PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')


def site_pdb():
    '''synthetic chain PDB with ligand LIG 300 + water HOH 501.

    The file lives in a temp dir that is kept (not deleted) until the
    process exits - the structure is re-read from disk by the pipeline.
    '''
    tmp_dir = tempfile.mkdtemp(prefix='sp_ligand_')
    return synthesize_site_pdb(Path(tmp_dir) / 'site.pdb')


def run_params(outdir, **kw):
    params = {
        'table_name': 'chain_ligand', 'host': 'test-host', 'verbose': 0,
        'overwrite': True,
        'output_pdb_dir': Path(outdir) if outdir else outdir,
        'search_cutoff': None, 'site_residues': None,
    }
    params.update(kw)
    return params


def site_row(pdb_file):
    return {'pdb_id': 'XXXX', 'chain_id': 'XXXX_A', 'chain_no': 'A',
            'relative_pdb': 'XXXX_A.pdb', 'pdb_file': pdb_file}


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestLigandNeighbors(TestCase):
    '''ligand_neighbors() wraps Neighbors.run() (no DB required)'''

    def test_ligand_neighbors_with_ligand(self):
        from pipelines.dynamics.chain_ligand import ligand_neighbors
        info = ligand_neighbors(site_pdb(), chain_no='A')
        self.assertEqual(info['n_sites'], 1)
        self.assertEqual(info['ligand_names'], 'LIG')
        self.assertEqual(info['n_neighbor_residues'], 4)
        self.assertEqual(info['n_fragments'], 2)
        self.assertEqual(len(info['neighbors']), 4)

    def test_ligand_neighbors_without_ligand(self):
        from pipelines.dynamics.chain_ligand import ligand_neighbors
        info = ligand_neighbors(PDB_2OX8, chain_no='A')
        self.assertEqual(info['n_sites'], 0)
        self.assertEqual(info['n_neighbor_residues'], 0)

    def test_ligand_neighbors_missing_file(self):
        from pipelines.dynamics.chain_ligand import ligand_neighbors
        self.assertIsNone(
            ligand_neighbors(os.path.join(DATA_DIR, 'nope.pdb')))

    def test_ligand_neighbors_chain_scoped_sites(self):
        from pipelines.dynamics.chain_ligand import ligand_neighbors
        # LIG300 / HOH501 sit in chain A: another chain has no sites
        info = ligand_neighbors(site_pdb(), chain_no='B')
        self.assertEqual(info['n_sites'], 0)
        self.assertEqual(info['n_neighbor_residues'], 0)

    def test_ligand_neighbors_model_no(self):
        from pipelines.dynamics.chain_ligand import ligand_neighbors
        # model_no=0 selects the (only) model of the file
        info = ligand_neighbors(site_pdb(), chain_no='A', model_no=0)
        self.assertEqual(info['n_sites'], 1)
        self.assertEqual(info['ligand_names'], 'LIG')


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestLigandPipelineRun(TestCase):
    '''run() exports contacts and yields chain_ligand records'''

    def test_run_yields_table_columns(self):
        from pipelines.dynamics.chain_ligand import run
        with tempfile.TemporaryDirectory() as tmp_dir:
            records = list(run([site_row(site_pdb())],
                run_params(tmp_dir), defaultdict(int)))
            self.assertEqual(len(records), 1)
            rec = dict(zip(TABLE_COLUMNS['chain_ligand'], records[0]))
            self.assertEqual(rec['pdb_id'], 'XXXX')
            self.assertEqual(rec['chain_id'], 'XXXX_A')
            self.assertEqual(rec['n_sites'], 1)
            self.assertEqual(rec['n_neighbor_residues'], 4)
            self.assertEqual(rec['n_fragments'], 2)
            self.assertEqual(rec['ligand_names'], 'LIG')
            # exports are written below output_pdb_dir and recorded
            for col in ('neighbor_csv', 'neighbor_pkl'):
                self.assertTrue((Path(tmp_dir) / rec[col]).is_file())
            contacts = pd.read_csv(Path(tmp_dir) / rec['neighbor_csv'])
            self.assertEqual(sorted(contacts['resseq'].tolist()),
                             [2, 3, 102, 103])
            with open(Path(tmp_dir) / rec['neighbor_pkl'], 'rb') as f:
                data = pickle.load(f)
            self.assertEqual(data['summary']['n_sites'], 1)

    def test_run_activity_site(self):
        from pipelines.dynamics.chain_ligand import run
        with tempfile.TemporaryDirectory() as tmp_dir:
            records = list(run([site_row(site_pdb())],
                run_params(tmp_dir, site_residues='A/2'),
                defaultdict(int)))
            rec = dict(zip(TABLE_COLUMNS['chain_ligand'], records[0]))
            self.assertEqual(rec['n_sites'], 1)
            self.assertEqual(rec['n_neighbor_residues'], 2)
            self.assertEqual(rec['n_fragments'], 2)

    def test_run_chain_without_site(self):
        from pipelines.dynamics.chain_ligand import run
        with tempfile.TemporaryDirectory() as tmp_dir:
            records = list(run([site_row(PDB_2OX8)],
                run_params(tmp_dir), defaultdict(int)))
            self.assertEqual(len(records), 1)
            rec = dict(zip(TABLE_COLUMNS['chain_ligand'], records[0]))
            self.assertEqual(rec['n_sites'], 0)
            self.assertIsNone(rec['neighbor_csv'])
            self.assertIsNone(rec['neighbor_pkl'])
            self.assertEqual(list(Path(tmp_dir).iterdir()), [])

    def test_run_skips_unparsable_files(self):
        from pipelines.dynamics.chain_ligand import run
        rows = [site_row(os.path.join(DATA_DIR, 'nope.pdb'))]
        self.assertEqual(
            list(run(rows, run_params(None), defaultdict(int))), [])


class TestLigandPipelineWiring(TestCase):
    '''table / command registration (no DB required)'''

    def test_table_columns_cover_chain_ligand(self):
        self.assertEqual(TABLE_COLUMNS['chain_ligand'],
                         ['pdb_id', 'chain_id', 'n_sites', 'ligand_names',
                          'n_neighbor_residues', 'n_fragments',
                          'min_dist', 'mean_dist', 'neighbor_residues',
                          'neighbor_csv', 'neighbor_pkl'])

    def test_main_pool_registers_chain_ligand(self):
        from src.app.main import pool
        self.assertEqual(pool['chain_ligand'], 'dynamics.chain_ligand')

    def test_sql_schema_defines_chain_ligand(self):
        sql_file = os.path.join(profiling_dir, 'src', 'sql', 'ligand.sql')
        with open(sql_file) as f:
            sql = f.read()
        self.assertIn('CREATE TABLE chain_ligand', sql)
        for col in ('n_sites', 'ligand_names', 'n_neighbor_residues',
                    'n_fragments', 'min_dist', 'mean_dist',
                    'neighbor_residues', 'neighbor_csv', 'neighbor_pkl'):
            self.assertIn(col, sql)


class TestParseSiteResidues(TestCase):
    '''site_residues CLI param parsing (no DB required)'''

    def test_parse_site_residues(self):
        from pipelines.dynamics.chain_ligand import parse_site_residues
        self.assertIsNone(parse_site_residues(None))
        self.assertIsNone(parse_site_residues(''))
        self.assertEqual(parse_site_residues('A/57,A/102'),
                         [('A', 57), ('A', 102)])
        self.assertEqual(parse_site_residues('B/12'), [('B', 12)])