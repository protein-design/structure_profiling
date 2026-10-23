'''
tests for src.helper.validate_data (shared row-validation helper of the
dynamics / pdb pipelines: chain_fixed, func_fixed, chain_relaxed,
func_relaxed, func_combo) with the database access mocked out:
- rows whose files exist are kept, missing ones are collected and deleted
- the deletion is chunked (chunk_size = 100) and must cover the WHOLE
  delete pool: the tail chunk (n = 250 -> chunks 100 + 100 + 50,
  n = 101 -> chunks 100 + 1) must be deleted too
  (regression: the old loop bounded range() at n - chunk_size, so the
  last <= chunk_size rows were never deleted)
- meta['delete_rows'] counts only successful delete_rows() calls
- nothing is deleted when every file exists
- pytest -s tests/simulation/test_validate_data.py
'''
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest import mock, skipUnless

project_dir = Path(__file__).resolve().parents[2]
HAS_ENV = (project_dir / '.env').is_file()


def make_rows(col_name, names):
    '''rows as returned by QueryComplex.list_data (pdb_id + col_name).'''
    return [{'pdb_id': name.replace('/', '_'), col_name: name}
            for name in names]


@skipUnless(HAS_ENV, 'requires local .env (output_dir etc.)')
class TestValidateData(unittest.TestCase):
    '''the chunked deletion of validate_data must cover the whole pool'''

    COL = 'fixed_pdb'
    TABLE = 'chain_fixed'

    def run_validate(self, names, existing=(), delete_result=True):
        '''create `existing` files in a tmp dir, run validate_data on
        `names` with mocked DB classes; returns (delete calls, meta).'''
        with tempfile.TemporaryDirectory() as tmp:
            for name in existing:
                path = Path(tmp) / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('ATOM\n')
            rows = make_rows(self.COL, names)
            calls = []

            def fake_delete(table_name, col_name, values):
                calls.append(list(values))
                return delete_result

            params = {'pdb_group': 0, 'table_name': self.TABLE,
                      'verbose': False}
            meta = defaultdict(int)
            # src.helper imports the DB classes into its own namespace
            from src.helper import validate_data
            with mock.patch('src.helper.QueryComplex') as qc_cls, \
                    mock.patch('src.helper.DeleteComplex') as dc_cls:
                qc_cls.return_value.list_data.return_value = rows
                dc_cls.return_value.delete_rows.side_effect = fake_delete
                validate_data(self.COL, params, meta, outdir=Path(tmp))
        return calls, meta

    def test_chunked_deletion_covers_tail(self):
        '''n = 250 -> chunks of 100 + 100 + 50, every value deleted once'''
        names = [f'missing_{i}.pdb' for i in range(250)]
        calls, meta = self.run_validate(names)
        self.assertEqual([len(c) for c in calls], [100, 100, 50])
        self.assertEqual([v for c in calls for v in c], names)
        self.assertEqual(meta['delete_rows'], 250)

    def test_deletion_just_over_one_chunk(self):
        '''n = 101 -> the single tail row must be deleted too'''
        names = [f'missing_{i}.pdb' for i in range(101)]
        calls, meta = self.run_validate(names)
        self.assertEqual([len(c) for c in calls], [100, 1])
        self.assertEqual(meta['delete_rows'], 101)

    def test_deletion_exact_one_chunk(self):
        '''n = 100 -> one full chunk'''
        names = [f'missing_{i}.pdb' for i in range(100)]
        calls, meta = self.run_validate(names)
        self.assertEqual([len(c) for c in calls], [100])
        self.assertEqual(meta['delete_rows'], 100)

    def test_deletion_below_one_chunk(self):
        '''n < chunk_size -> a single partial chunk'''
        names = [f'missing_{i}.pdb' for i in range(99)]
        calls, meta = self.run_validate(names)
        self.assertEqual([len(c) for c in calls], [99])
        self.assertEqual(meta['delete_rows'], 99)

    def test_existing_files_are_kept(self):
        '''only rows whose file is missing are deleted'''
        names = [f'row_{i}.pdb' for i in range(250)]
        existing = [f'row_{i}.pdb' for i in range(200)]
        calls, meta = self.run_validate(names, existing=existing)
        deleted = [v for c in calls for v in c]
        self.assertEqual(deleted, [f'row_{i}.pdb' for i in range(200, 250)])
        self.assertEqual(meta['delete_rows'], 50)

    def test_nothing_missing_no_deletion(self):
        '''every file exists -> delete_rows is never called'''
        names = [f'row_{i}.pdb' for i in range(120)]
        calls, meta = self.run_validate(names, existing=names)
        self.assertEqual(calls, [])
        self.assertNotIn('delete_rows', meta)

    def test_failed_delete_not_counted(self):
        '''a falsy delete_rows() result must not update the meta counter'''
        names = [f'missing_{i}.pdb' for i in range(150)]
        calls, meta = self.run_validate(names, delete_result=False)
        self.assertEqual([len(c) for c in calls], [100, 50])
        self.assertNotIn('delete_rows', meta)


if __name__ == '__main__':
    unittest.main()
