'''
tests for pipelines/test_env.py: read-only environment diagnostics.
No database required: the mysql probe is mocked out (CI hosts have no
COMPLEX2), everything else only reads the local .env loaded by tests/helper.py.
- pytest -s tests/app/test_env_pipeline.py
'''
import tempfile
from collections import defaultdict
from pathlib import Path
from unittest import TestCase, mock

from pipelines.test_env import (
    MISS,
    OK,
    OPTIONAL_MODULES,
    REQUIRED_MODULES,
    SECRET_MARKS,
    WARN,
    check_binaries,
    check_cuda,
    check_env_file,
    check_imports,
    run,
)

from ..helper import *


def statuses(rows):
    '''map report rows to {item: status}'''
    return {row['item']: row['status'] for row in rows}


class TestCheckImports(TestCase):

    def test_python_and_required_modules_ok(self):
        '''hard dependencies of pyproject.toml must import'''
        found = statuses(check_imports())
        self.assertIn('python', found)
        for name in REQUIRED_MODULES:
            self.assertEqual(found[name], OK, f'{name} is a hard dependency')

    def test_optional_modules_never_miss(self):
        '''optional extras (torch etc.) may be absent but are never MISS'''
        found = statuses(check_imports())
        for name in OPTIONAL_MODULES:
            self.assertIn(found[name], (OK, WARN))


class TestCheckEnvFile(TestCase):
    '''running the test suite already requires a local .env'''

    def setUp(self):
        self.rows = check_env_file()
        self.found = statuses(self.rows)

    def test_env_path_resolved(self):
        self.assertEqual(self.found['.env path'], OK)

    def test_required_keys_present(self):
        self.assertEqual(self.found['required keys'], OK)

    def test_output_dir_writable(self):
        self.assertEqual(self.found['output_dir'], OK)

    def test_secrets_masked(self):
        '''the report must never print a row for a secret key
        (masking is by key: a secret value may coincide with a shown value)'''
        for row in self.rows:
            for mark in SECRET_MARKS:
                self.assertNotIn(mark, row['item'].upper())


class TestCheckCuda(TestCase):

    def test_at_least_one_platform_listed(self):
        '''openmm always registers the CPU platform'''
        platform_rows = [r for r in check_cuda() if r['item'].startswith('openmm platform')]
        self.assertTrue(platform_rows)
        self.assertTrue(all(r['status'] == OK for r in platform_rows))

    def test_cuda_rows_are_advisory(self):
        '''cuda/torch availability is machine dependent: only OK or WARN'''
        found = statuses(check_cuda())
        self.assertIn(found['openmm CUDA'], (OK, WARN))
        self.assertIn(found['pytorch cuda'], (OK, WARN))


class TestCheckBinaries(TestCase):

    def test_mkdssp_reported(self):
        rows = check_binaries()
        self.assertEqual(rows[0]['item'], 'mkdssp')
        self.assertIn(rows[0]['status'], (OK, WARN))

    def _ccd_rows(self, tmp, ccd_file=None):
        '''check_binaries() against a hermetic fake data dir and binary'''
        if ccd_file is not None:
            Path(tmp, ccd_file).write_text('fake ccd', encoding='utf-8')
        with mock.patch('pipelines.test_env.CCD_DIRS', (tmp,)), \
                mock.patch('pipelines.test_env.shutil.which',
                           return_value='/usr/bin/mkdssp'):
            return statuses(check_binaries())

    def test_ccd_gz_found_is_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            found = self._ccd_rows(tmp, 'components.cif.gz')
        self.assertEqual(found['mkdssp CCD'], OK)

    def test_ccd_plain_found_is_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            found = self._ccd_rows(tmp, 'components.cif')
        self.assertEqual(found['mkdssp CCD'], OK)

    def test_ccd_missing_is_warn(self):
        '''absent CCD only degrades to WARN, it is not a hard MISS'''
        with tempfile.TemporaryDirectory() as tmp:
            found = self._ccd_rows(tmp)
        self.assertEqual(found['mkdssp CCD'], WARN)


class TestRun(TestCase):

    def test_run_tallies_every_row(self):
        '''run() counts OK/WARN/MISS for all printed rows (mysql mocked)'''
        with mock.patch('pipelines.test_env.check_mysql',
                        return_value=[{'item': 'mysql', 'status': OK, 'detail': 'mocked'}]):
            meta = defaultdict(int)
            sections = run({}, meta)
        n_rows = sum(len(rows) for _, rows in sections)
        self.assertEqual(sum(meta[s] for s in (OK, WARN, MISS)), n_rows)
        self.assertGreater(meta[OK], 0)
        self.assertNotIn('FAIL', meta)