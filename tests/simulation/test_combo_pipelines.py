'''
tests for the combo split pipelines moved from bio_omics:
combo2..5 and func_combo (shared), plus the ported helpers
QueryComplex.pdb_chains, Pdb.to_chain_pdb and ProcessSeq.
Importing the pipeline modules must not require a DB connection,
only a local .env (output_dir etc.) used by src.helper at import time.
'''
import os
import tempfile
import unittest
from pathlib import Path
from unittest import skipUnless

project_dir = Path(__file__).resolve().parents[2]
HAS_ENV = (project_dir / '.env').is_file()
PDB_FILE = project_dir / 'tests' / 'data' / 'PDB' / '2OX8_A.pdb'


class TestProcessSeq(unittest.TestCase):
    '''ported colabfold exporters (no DB required)'''

    def test_to_colabfold_csv(self):
        from structure_profiling import ProcessSeq
        with tempfile.TemporaryDirectory() as tmp_dir:
            outfile = os.path.join(tmp_dir, 'colabfold.csv')
            data = {'complex1': ['AAA', 'DDD']}
            assert ProcessSeq.to_colabfold_csv(data, outfile) == 'succeed'
            with open(outfile) as f:
                assert f.read() == 'id,sequence\ncomplex1,AAA:DDD\n'

    def test_to_colabfold_fasta(self):
        from structure_profiling import ProcessSeq
        with tempfile.TemporaryDirectory() as tmp_dir:
            outfile = os.path.join(tmp_dir, 'colabfold.faa')
            data = {'complex1': ['AAA', 'DDD']}
            assert ProcessSeq.to_colabfold_fasta(data, outfile) == 'succeed'
            with open(outfile) as f:
                assert f.read() == '>complex1\nAAA:\nDDD'


@skipUnless(PDB_FILE.is_file(), 'requires tests/data/PDB/2OX8_A.pdb')
class TestPdbToChainPdb(unittest.TestCase):
    '''ported Pdb.to_chain_pdb (no DB required)'''

    def test_to_chain_pdb(self):
        from structure_profiling import Pdb
        p = Pdb(PDB_FILE, None, False)
        self.assertTrue(p.structure)
        with tempfile.TemporaryDirectory() as tmp_dir:
            outfile = Path(tmp_dir) / 'combo.pdb'
            status = p.to_chain_pdb(0, ['A'], outfile)
            self.assertTrue(status)
            self.assertTrue(outfile.is_file())
            with open(outfile) as f:
                lines = [line for line in f if line.startswith(('ATOM', 'HETATM'))]
            self.assertTrue(len(lines) > 0)


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestComboPipelineImports(unittest.TestCase):

    def test_func_combo(self):
        from src.pipelines.pdb import func_combo
        self.assertTrue(callable(func_combo.retrieve_data))
        self.assertTrue(callable(func_combo.build_data))

    def test_combo_pipelines(self):
        from src.pipelines.pdb import combo2, combo3, combo4, combo5
        for module in (combo2, combo3, combo4, combo5):
            self.assertTrue(callable(module.retrieve_data))
            self.assertTrue(callable(module.build_data))

    def test_query_complex_has_pdb_chains(self):
        from src.structure_profiling.database.query_complex import QueryComplex
        self.assertTrue(callable(QueryComplex.pdb_chains))

    def test_table_columns_cover_combo_tables(self):
        from src.structure_profiling.database.build_complex import TABLE_COLUMNS
        for table in ('combo2', 'combo3', 'combo4', 'combo5'):
            self.assertEqual(TABLE_COLUMNS[table][-1], 'combo_fasta')