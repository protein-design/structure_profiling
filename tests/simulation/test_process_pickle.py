'''
tests for common/process_pickle.py (ProcessPickle):
- save_pickle / load_pickle round trip and the 'succeed' / 'fail'
  statuses (the unreachable duplicated 'return fail' was removed)
- save_data merges into an existing pickle
- load_pickle returns the default for missing / corrupted files
- is_good() is False for empty, missing or corrupted pickles
- pytest -s tests/simulation/test_process_pickle.py
'''
import tempfile
import unittest
from pathlib import Path

from structure_profiling import ProcessPickle


class TestProcessPickle(unittest.TestCase):

    def test_save_load_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            pfile = Path(tmp) / 'data.pkl'
            data = {'chain1': [1.0, 2.0], 'chain2': {'a': 1}}
            self.assertEqual(ProcessPickle(pfile).save_pickle(data),
                             'succeed')
            self.assertEqual(ProcessPickle(pfile).load_pickle(), data)

    def test_save_data_merges_existing(self):
        with tempfile.TemporaryDirectory() as tmp:
            pfile = Path(tmp) / 'data.pkl'
            pp = ProcessPickle(pfile)
            self.assertEqual(pp.save_data({'a': 1}), 'succeed')
            self.assertEqual(pp.save_data({'b': 2}), 'succeed')
            self.assertEqual(pp.load_pickle(), {'a': 1, 'b': 2})

    def test_load_missing_returns_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            pfile = Path(tmp) / 'absent.pkl'
            pp = ProcessPickle(pfile)
            self.assertIsNone(pp.load_pickle())
            self.assertEqual(pp.load_pickle(default={}), {})
            self.assertFalse(pp.is_good())

    def test_save_failure_returns_fail(self):
        '''an unwritable target path must yield 'fail', not raise'''
        with tempfile.TemporaryDirectory() as tmp:
            blocker = Path(tmp) / 'blocker'
            blocker.write_text('a plain file, not a directory')
            pfile = blocker / 'sub' / 'data.pkl'
            pp = ProcessPickle(pfile, verbose=True)
            self.assertEqual(pp.save_pickle({'a': 1}), 'fail')
            self.assertFalse(pp.is_good())

    def test_corrupted_pickle_returns_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            pfile = Path(tmp) / 'corrupt.pkl'
            pfile.write_bytes(b'this is not a pickle')
            pp = ProcessPickle(pfile)
            self.assertIsNone(pp.load_pickle())
            self.assertEqual(pp.load_pickle(default='fallback'), 'fallback')
            self.assertFalse(pp.is_good())

    def test_is_good_after_save(self):
        with tempfile.TemporaryDirectory() as tmp:
            pfile = Path(tmp) / 'data.pkl'
            pp = ProcessPickle(pfile)
            pp.save_pickle({'a': 1})
            self.assertTrue(pp.is_good())

    def test_accepts_str_path(self):
        '''str inputs are coerced to Path (see __init__ comment)'''
        with tempfile.TemporaryDirectory() as tmp:
            pfile = str(Path(tmp) / 'data.pkl')
            self.assertEqual(ProcessPickle(pfile).save_pickle({'a': 1}),
                             'succeed')
            self.assertEqual(ProcessPickle(pfile).load_pickle(), {'a': 1})


if __name__ == '__main__':
    unittest.main()
