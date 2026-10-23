'''
tests for structure_profiling/common/clean.py (Clean.clean_msa):
the af_msa pipeline drives colabfold_search, which spawns the mmseqs
workflow (colabfold_search -> blastpgp.sh -> mmseqs). Interrupted runs
leave orphaned mmseqs processes behind (holding GPU/RAM): clean_msa()
must kill those, while the mmseqs processes of running MSA jobs (a live
colabfold_search driver above them) and unrelated mmseqs processes
must survive. NB: the real colabfold_search driver runs as
'python3.13 .../envs/.../colabfold_search ...' (interpreter as argv0),
which the driver detection must recognize.
- pytest -s tests/simulation/test_clean.py
'''
import os
import shutil
import subprocess
import sys
import tempfile
import time
from unittest import TestCase

from src.structure_profiling import Clean

PY = sys.executable
MARKER = 'colabfold_search_99_TEST_0_X-Y_1_2'
CODE = f'import time; time.sleep(60) # {MARKER}'


class TestCleanMsa(TestCase):

    def setUp(self):
        # a fake colabfold_search driver script: the process is spawned
        # as 'python3.13 .../colabfold_search ...' exactly like the real
        # conda entrypoint (interpreter as argv0, script path as argv1)
        tmp_dir = tempfile.mkdtemp(prefix='sp_clean_test_')
        self.addCleanup(shutil.rmtree, tmp_dir, ignore_errors=True)
        driver_path = os.path.join(tmp_dir, 'colabfold_search')
        with open(driver_path, 'w') as fh:
            fh.write(
                "import subprocess, time;"
                f"subprocess.Popen(['mmseqs', '-c', {CODE!r}],"
                f"executable={PY!r});"
                "time.sleep(60)"
            )
        os.chmod(driver_path, 0o755)

        # orphaned mmseqs of an interrupted run: MSA-marked command line,
        # no live colabfold_search driver above it
        self.orphan = subprocess.Popen(
            ['mmseqs', '-c', CODE], executable=PY)
        # running MSA job: a live colabfold_search driver (real process
        # signature) whose child is a MSA-marked mmseqs process
        self.driver = subprocess.Popen([PY, driver_path, '--threads', '8'])
        # workflow arguments referencing the colabfold_search tmpdir
        # (like .../colabfold_search_36_<input>/qdb): their basenames
        # start with 'colabfold_search' but must never count as a driver,
        # or orphaned subtrees would protect themselves
        self.workflow = subprocess.Popen(
            ['mmseqs', '-c', CODE, os.path.join(tmp_dir, MARKER)],
            executable=PY)
        # an unrelated mmseqs run (no colabfold markers): untouched
        self.foreign = subprocess.Popen(
            ['mmseqs', '-c', 'import time; time.sleep(60)'], executable=PY)
        # let the process trees settle
        time.sleep(1)

    def tearDown(self):
        for proc in (self.orphan, self.driver, self.workflow, self.foreign):
            if proc.poll() is None:
                proc.kill()
        # reap the mmseqs children orphaned by the kills above
        Clean.clean_msa(grace=1.0)

    def _alive(self, pid: int) -> bool:
        try:
            with open(f'/proc/{pid}/stat') as fh:
                return fh.read().rsplit(')', 1)[1].split()[0] != 'Z'
        except (FileNotFoundError, ProcessLookupError):
            return False

    def test_kills_orphans_keeps_running_msa_and_foreign(self):
        self.assertIsNone(self.orphan.poll())
        self.assertIsNone(self.driver.poll())
        self.assertIsNone(self.foreign.poll())

        killed = Clean.clean_msa(grace=1.0, verbose=True)

        # the orphaned mmseqs of the interrupted run is killed
        self.assertIn(self.orphan.pid, killed)
        self.assertFalse(self._alive(self.orphan.pid))
        # the running MSA job survives: driver and its mmseqs child
        self.assertIsNone(self.driver.poll())
        children = subprocess.run(
            ['pgrep', '-P', str(self.driver.pid)],
            check=False, capture_output=True, text=True).stdout.split()
        self.assertEqual(1, len(children))
        self.assertTrue(self._alive(int(children[0])))
        # the driverless workflow tree (tmpdir-named arguments) is killed
        self.assertIn(self.workflow.pid, killed)
        self.assertFalse(self._alive(self.workflow.pid))
        # the unrelated mmseqs run is untouched
        self.assertIsNone(self.foreign.poll())
        self.assertTrue(self._alive(self.foreign.pid))
