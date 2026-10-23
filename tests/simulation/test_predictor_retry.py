'''
tests for structure_profiling/predict/predictor.py (colabfold_search retry):
the GPU ungapped prefilter needs a large workspace (cudaMalloc ~6 GiB for
uniref30) and dies with std::bad_alloc when video memory is held by other
processes. The retry must (a) wait for free VRAM and retry on GPU,
(b) fall back to a CPU search passing --prefilter-mode 1 explicitly,
because the local *.idx databases only carry ungapped-prefilter data and
colabfold's CPU default (prefilter-mode 0, kmer) always fails with
"Index was not built with prefilter support".
- pytest -s tests/simulation/test_predictor_retry.py
'''
from unittest import TestCase, mock

from src.structure_profiling import Predictor

ROW = {'query_faa': '/tmp/q.fasta', 'outdir': '/tmp/out'}
PARAMS = {
    'gpu': True,
    'colabfold_db': '/db',
    'colabfold_bin_dir': '/bin',
    'tools_dir': '/tools',
    'gpu_vram_wait': 0.0,
}

def _cmd_str(call) -> str:
    return ' '.join(str(part) for part in call.args[1])


class TestColabfoldSearchRetry(TestCase):

    def run_search(self, statuses, free_mib_values, params=None):
        '''mock colabfold_subprocess and gpu_free_mib, run the search'''
        params = dict(PARAMS if params is None else params)
        with mock.patch.object(
            Predictor, 'colabfold_subprocess', side_effect=statuses,
        ) as proc, mock.patch.object(
            Predictor, 'gpu_free_mib', side_effect=free_mib_values,
        ), mock.patch(
            'src.structure_profiling.predict.predictor.time.sleep',
        ) as sleep:
            status = Predictor.colabfold_search(ROW, params)
        return status, proc.call_args_list, sleep.call_args_list

    def test_gpu_success_needs_no_retry(self):
        status, calls, sleeps = self.run_search([0], [9000])
        self.assertEqual(0, status)
        self.assertEqual(1, len(calls))
        self.assertIn('--gpu 1', _cmd_str(calls[0]))
        self.assertEqual([], sleeps)

    def test_waits_for_vram_then_retries_on_gpu(self):
        # first attempt dies; VRAM is busy twice before it frees up
        params = dict(PARAMS, gpu_vram_wait=30)
        status, calls, sleeps = self.run_search(
            [1, 0], [1000, 1000, 9000, 9000], params,
        )
        self.assertEqual(0, status)
        self.assertEqual(2, len(calls))
        self.assertEqual(['15', '15'], [str(c.args[0]) for c in sleeps])
        self.assertIn('--gpu 1', _cmd_str(calls[0]))
        self.assertIn('--gpu 1', _cmd_str(calls[1]))

    def test_cpu_fallback_uses_prefilter_mode_1(self):
        # VRAM frees up but the GPU retry still dies: CPU search must
        # pass --gpu 0 --prefilter-mode 1 (ungapped prefilter)
        status, calls, sleeps = self.run_search(
            [1, 1, 7], [9000, 9000],
        )
        self.assertEqual(7, status)
        self.assertEqual(3, len(calls))
        self.assertEqual([], sleeps)
        cpu_cmd = _cmd_str(calls[2])
        self.assertIn('--gpu 0', cpu_cmd)
        self.assertIn('--prefilter-mode 1', cpu_cmd)
        self.assertNotIn('--gpu 1', cpu_cmd)
        self.assertIn('/tmp/q.fasta', cpu_cmd)

    def test_cpu_fallback_when_vram_never_frees(self):
        # nvidia-smi reports almost no free VRAM the whole time: after the
        # wait timeout the CPU search is the only chance left
        params = dict(PARAMS, gpu_vram_wait=15)
        status, calls, sleeps = self.run_search(
            [1, 3], [1000, 1000, 1000], params,
        )
        self.assertEqual(3, status)
        self.assertEqual(2, len(calls))
        self.assertEqual(['15'], [str(c.args[0]) for c in sleeps])
        self.assertIn('--prefilter-mode 1', _cmd_str(calls[1]))

    def test_gpu_free_mib_parses_nvidia_smi(self):
        fake = mock.Mock(returncode=0, stdout='9892\n', stderr='')
        with mock.patch(
            'src.structure_profiling.predict.predictor.subprocess.run',
            return_value=fake,
        ):
            self.assertEqual(9892, Predictor.gpu_free_mib())
        fake = mock.Mock(returncode=1, stdout='', stderr='driver down')
        with mock.patch(
            'src.structure_profiling.predict.predictor.subprocess.run',
            return_value=fake,
        ):
            self.assertEqual(0, Predictor.gpu_free_mib())

    def test_no_retry_without_gpu_param(self):
        status, calls, _ = self.run_search(
            [1], [], dict(PARAMS, gpu=False),
        )
        self.assertEqual(1, status)
        self.assertEqual(1, len(calls))
        self.assertNotIn('--gpu', _cmd_str(calls[0]))