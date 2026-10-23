import builtins
import signal
import functools

# subclass the built-in TimeoutError so that a caller can catch a timeout
# with either this class or the built-in one ('except TimeoutError' in
# openmm_pdbfixer.py / openmm_simulation.py must keep working)
class TimeoutError(builtins.TimeoutError):
    pass

def timeout(minutes):
    seconds = int(minutes * 60)
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            def timeout_handler(signum, frame):
                raise TimeoutError(f"{func.__name__} timed out after {minutes} minutes")
            
            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            return result
        return wrapper
    return decorator