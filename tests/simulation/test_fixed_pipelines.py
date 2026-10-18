'''
smoke tests for the fixed pipelines moved from bio_omics:
chain_fixed, func_fixed (shared) and combo2..5_fixed.
Importing them must not require a DB connection, only a local .env
(output_dir etc.) used by src.helper at import time.
'''
import unittest
from pathlib import Path
from unittest import skipUnless

project_dir = Path(__file__).resolve().parents[2]
HAS_ENV = (project_dir / '.env').is_file()


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestFixedPipelineImports(unittest.TestCase):

    def test_chain_fixed(self):
        from src.pipelines.dynamics import chain_fixed
        self.assertTrue(callable(chain_fixed.retrieve_pdb))
        self.assertTrue(callable(chain_fixed.fix))

    def test_func_fixed(self):
        from src.pipelines.dynamics import func_fixed
        self.assertTrue(callable(func_fixed.retrieve_pdb))
        self.assertTrue(callable(func_fixed.fix))

    def test_combo_fixed(self):
        from src.pipelines.dynamics import (combo2_fixed, combo3_fixed,
            combo4_fixed, combo5_fixed)
        for module in (combo2_fixed, combo3_fixed, combo4_fixed, combo5_fixed):
            self.assertTrue(callable(module.retrieve_pdb))
            self.assertTrue(callable(module.fix))

    def test_query_complex_has_retrieve_combo(self):
        from src.structure_profiling.database.query_complex import QueryComplex
        self.assertTrue(callable(QueryComplex.retrieve_combo))

    def test_table_columns_cover_fixed_tables(self):
        from src.structure_profiling.database.build_complex import TABLE_COLUMNS
        for table in ('chain_fixed', 'combo2_fixed', 'combo3_fixed',
                      'combo4_fixed', 'combo5_fixed'):
            self.assertEqual(TABLE_COLUMNS[table][-1], 'fixed_pdb')
