'''
tests for pipelines/predict/func_summary.py with the MMalign binary
mocked out (shared by af_rmsd.py and af_tmalign.py):
- calculate_tmalign yields one record per aligned structure pair and
  counts it under f'{pdb_status}_{category}_succeed'
- a failed alignment (TMalign returns the empty dict {}) is counted
  under f'{pdb_status}_{category}_fail' and yields NO record
  (regression: the old code counted into the bogus key
  f'{[prefix]}_fail' - a list interpolated into the string - and then
  still yielded a NULL-score record)
- pairs with a missing pdb file are skipped without counting
- calculate_rmsd pairs native/model paths and keeps only existing files
  (EvalPred mocked, no superposition actually computed)
- pytest -s tests/predict/test_func_summary.py
'''
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest import mock, skipUnless

project_dir = Path(__file__).resolve().parents[2]
HAS_ENV = (project_dir / '.env').is_file()


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestCalculateTmalign(unittest.TestCase):
    '''TM-align record builder: success / failure / missing-file paths'''

    def make_rows(self, tmp, native_exists=True, model_exists=True):
        '''one input row with (optionally) existing pdb files'''
        native_dir = Path(tmp) / 'native'
        pred_dir = Path(tmp) / 'pred'
        native_dir.mkdir()
        pred_dir.mkdir()
        native_pdb = native_dir / 'chain.pdb'
        model_pdb = pred_dir / 'model.pdb'
        if native_exists:
            native_pdb.write_text('ATOM\n')
        if model_exists:
            model_pdb.write_text('ATOM\n')
        return [{
            'input_id': 'chain1',
            'pdb_status': 'raw',
            'native_outdir': str(native_dir),
            'native_pdb': 'chain.pdb',
            'pred_outdir': str(pred_dir),
            'pred_pdb_name': 'model.pdb',
        }]

    def params(self, tmp):
        return {'pdb_status': 'raw', 'category': 'chain',
                'mmalign_dir': Path(tmp)}

    def test_success_yields_record_and_counts_succeed(self):
        from src.pipelines.predict.func_summary import calculate_tmalign
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.make_rows(tmp)
            scores = {'RMSD': '1.53', 'TM1': '0.92', 'TM2': '0.90'}
            with mock.patch(
                    'src.pipelines.predict.func_summary.EvalPred') as ep_cls:
                ep_cls.return_value.TMalign.return_value = scores
                meta = defaultdict(int)
                recs = list(calculate_tmalign(rows, self.params(tmp), meta))
            self.assertEqual(len(recs), 1)
            input_id, pdb_status, rmsd, tm1, tm2, native, model = recs[0]
            self.assertEqual(input_id, 'chain1')
            self.assertEqual(pdb_status, 'raw')
            self.assertEqual((rmsd, tm1, tm2), ('1.53', '0.92', '0.90'))
            self.assertEqual(native, 'chain.pdb')
            self.assertEqual(model, 'model.pdb')
            self.assertEqual(meta['raw_chain_succeed'], 1)
            self.assertEqual(meta['raw_chain_fail'], 0)
            # no bogus list-interpolated key may appear
            for key in meta:
                self.assertNotIn('[', key)

    def test_failure_counts_fail_and_yields_nothing(self):
        from src.pipelines.predict.func_summary import calculate_tmalign
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.make_rows(tmp)
            with mock.patch(
                    'src.pipelines.predict.func_summary.EvalPred') as ep_cls:
                # EvalPred.TMalign returns {} when the binary fails
                ep_cls.return_value.TMalign.return_value = {}
                meta = defaultdict(int)
                recs = list(calculate_tmalign(rows, self.params(tmp), meta))
            self.assertEqual(recs, [])
            self.assertEqual(meta['raw_chain_fail'], 1)
            self.assertEqual(meta['raw_chain_succeed'], 0)
            for key in meta:
                self.assertNotIn('[', key)

    def test_missing_native_pdb_skipped_without_counting(self):
        from src.pipelines.predict.func_summary import calculate_tmalign
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.make_rows(tmp, native_exists=False)
            with mock.patch(
                    'src.pipelines.predict.func_summary.EvalPred') as ep_cls:
                ep_cls.return_value.TMalign.return_value = {
                    'RMSD': '1.53', 'TM1': '0.92', 'TM2': '0.90'}
                meta = defaultdict(int)
                recs = list(calculate_tmalign(rows, self.params(tmp), meta))
            self.assertEqual(recs, [])
            self.assertNotIn('raw_chain_succeed', meta)
            self.assertNotIn('raw_chain_fail', meta)


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestCalculateRmsd(unittest.TestCase):
    '''RMSD record builder: path pairing and file filtering'''

    def test_pairs_kept_and_missing_filtered(self):
        from src.pipelines.predict.func_summary import calculate_rmsd
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'pdb').mkdir()
            (tmp / 'pred').mkdir()
            (tmp / 'pdb' / 'native.pdb').write_text('ATOM\n')
            (tmp / 'pred' / 'model.pdb').write_text('ATOM\n')
            good = {
                'input_id': 'chain1', 'pdb_status': 'raw',
                'native_outdir': str(tmp / 'pdb'),
                'native_pdb': 'native.pdb',
                'pred_outdir': str(tmp / 'pred'),
                'pred_pdb_name': 'model.pdb',
            }
            # native file of this row does not exist -> must be filtered
            broken = dict(good, input_id='chain2', native_pdb='nope.pdb')
            meta = defaultdict(int)
            with mock.patch(
                    'src.pipelines.predict.func_summary.EvalPred') as ep_cls:
                ep_cls.return_value.batch_rmsd.side_effect = \
                    lambda input_rows: iter([
                        dict(r, rmsd=0.75) for r in input_rows])
                recs = list(calculate_rmsd([good, broken], meta))
            self.assertEqual(meta['pred_pdb_files'], 1)
            self.assertEqual(len(recs), 1)
            input_id, pdb_status, rmsd, native, model = recs[0]
            self.assertEqual(input_id, 'chain1')
            self.assertEqual(pdb_status, 'raw')
            self.assertEqual(rmsd, 0.75)
            self.assertEqual(native, 'native.pdb')
            self.assertEqual(model, 'model.pdb')


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestModelPdbQueries(unittest.TestCase):
    '''query builders return (pdb_status, category, query) triples'''

    def base_params(self, tmp):
        return {'pdb_group': 0, 'params_table_name': 'af_params',
                'summary_table_name': 'af_summary',
                'table_name': 'af_rmsd', 'output_pdb_dir': Path(tmp),
                'simulate_dir': Path(tmp) / 'simulate'}

    def test_chain_query_statuses(self):
        from src.pipelines.predict.func_summary import (
            chain_raw_model_pdb, chain_fixed_model_pdb,
            chain_relaxed_model_pdb)
        with tempfile.TemporaryDirectory() as tmp:
            params = self.base_params(tmp)
            self.assertEqual(chain_raw_model_pdb(params)[0:2],
                             ('raw', 'chain'))
            self.assertEqual(chain_fixed_model_pdb(params)[0:2],
                             ('fixed', 'chain'))
            self.assertEqual(chain_relaxed_model_pdb(params)[0:2],
                             ('relaxed', 'chain'))

    def test_combo_query_status_and_category(self):
        from src.pipelines.predict.func_summary import combo_raw_model_pdb
        with tempfile.TemporaryDirectory() as tmp:
            params = self.base_params(tmp)
            pdb_status, category, query = combo_raw_model_pdb(params, 'combo3')
            self.assertEqual((pdb_status, category), ('raw', 'combo3'))
            self.assertIn('LEFT JOIN combo3', query)
            self.assertIn('c.combo_pdb', query)


if __name__ == '__main__':
    unittest.main()
