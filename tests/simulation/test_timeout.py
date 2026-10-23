'''
unit tests for the custom timeout decorator and its TimeoutError
- pytest -s tests/simulation/test_timeout.py
'''
import builtins
import signal
import time

from ..helper import *
from src.structure_profiling.structure.timeout import TimeoutError, timeout


@ddt
class TestTimeout(TestCase):

    def test_fast_function_passes(self):
        @timeout(2)
        def fast(a, b=0):
            'adds things'
            return a + b
        self.assertEqual(fast(1, b=2), 3)
        # functools.wraps keeps the identity of the decorated function
        self.assertEqual(fast.__name__, 'fast')
        self.assertEqual(fast.__doc__, 'adds things')

    def test_timeout_raises_custom_error(self):
        # the decorator's argument is MINUTES (0.05 min = 3 s)
        @timeout(0.05)
        def slow():
            time.sleep(4)
            return 'done'
        with self.assertRaises(TimeoutError):
            slow()

    def test_timeout_is_builtin_compatible(self):
        # openmm_pdbfixer / openmm_simulation rely on the custom
        # TimeoutError being catchable as the built-in one
        @timeout(0.05)
        def slow():
            time.sleep(4)
        with self.assertRaises(builtins.TimeoutError):
            slow()

    def test_timeout_message_names_function(self):
        @timeout(0.05)
        def slow_job():
            time.sleep(4)
        try:
            slow_job()
            self.fail('TimeoutError was not raised')
        except TimeoutError as err:
            self.assertIn('slow_job', str(err))
            self.assertIn('timed out', str(err))

    def test_no_pending_alarm_after_call(self):
        # the previous alarm state must be restored when the call succeeds
        @timeout(2)
        def quick():
            return 42
        self.assertEqual(quick(), 42)
        self.assertEqual(signal.alarm(0), 0)
