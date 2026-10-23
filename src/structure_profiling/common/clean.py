'''
Clean the leftover processes of interrupted pipeline runs.
The af_msa pipeline drives colabfold_search, which spawns the mmseqs
workflow (colabfold_search -> blastpgp.sh -> mmseqs ...). When a run is
killed, the colabfold_search driver dies but its mmseqs children survive
and keep holding GPU/RAM. clean_msa() kills those unused leftovers,
while the processes of running MSA jobs are kept.
'''
import os
import signal
import time


class Clean:

    # a process belongs to a colabfold MSA run when its own command line
    # references the colabfold_search working dirs or the mmseqs workflow
    MSA_MARKERS = ('colabfold_search', 'blastpgp.sh')
    # a live colabfold_search driver protects its whole subtree
    MSA_DRIVER = 'colabfold_search'

    @staticmethod
    def proc_stats() -> dict:
        '''
        Snapshot of the running processes:
        {pid: {'ppid': int, 'state': str, 'argv': [str, ...]}}
        Processes disappearing during the scan are skipped.
        '''
        stats = {}
        for entry in os.listdir('/proc'):
            if not entry.isdigit():
                continue
            pid = int(entry)
            try:
                with open(f'/proc/{pid}/stat') as fh:
                    stat = fh.read()
                # comm may contain spaces and brackets: parse after last ')'
                fields = stat[stat.rindex(')') + 2:].split()
                ppid, state = int(fields[1]), fields[0]
                with open(f'/proc/{pid}/cmdline', 'rb') as fh:
                    argv = [arg for arg
                            in fh.read().decode(errors='replace').split('\0')
                            if arg]
            except (FileNotFoundError, ProcessLookupError,
                    PermissionError, ValueError):
                # terminated during the scan or not accessible
                continue
            stats[pid] = {'ppid': ppid, 'state': state, 'argv': argv}
        return stats

    @staticmethod
    def _in_msa_tree(pid: int, stats: dict) -> bool:
        '''the process itself belongs to a colabfold MSA run'''
        seen = set()
        while pid > 0 and pid not in seen:
            seen.add(pid)
            proc = stats.get(pid)
            if proc is None:
                return False
            if any(marker in ' '.join(proc['argv'])
                    for marker in Clean.MSA_MARKERS):
                return True
            pid = proc['ppid']
        return False

    @staticmethod
    def _has_live_driver(pid: int, stats: dict) -> bool:
        '''a live colabfold_search driver among the strict ancestors'''
        proc = stats.get(pid)
        pid = proc['ppid'] if proc else 0
        seen = set()
        while pid > 1 and pid not in seen:
            seen.add(pid)
            proc = stats.get(pid)
            if proc is None:
                # parent is gone: orphaned process subtree
                return False
            # colabfold_search runs as 'python3.13 .../colabfold_search ...',
            # so any argv element with that basename marks the driver.
            # Exact basename match: a workflow tmpdir argument like
            # .../colabfold_search_36_<input>_... must not count as
            # a driver, or orphaned subtrees would protect themselves.
            if any(os.path.basename(arg) == Clean.MSA_DRIVER
                    for arg in proc['argv']):
                return True
            pid = proc['ppid']
        return False

    @staticmethod
    def clean_msa(verbose: bool = False, grace: float = 3.0) -> list:
        '''
        Kill the unused mmseqs processes left over by interrupted
        colabfold_search runs: mmseqs and its blastpgp.sh workflow without
        a live colabfold_search driver. The processes of running MSA jobs
        (subtree of a live colabfold_search) are kept, as well as mmseqs
        runs of other tools and other users.
        Returns the list of killed pids.
        '''
        stats = Clean.proc_stats()
        uid = os.getuid()
        targets = []
        kept = []
        for pid, proc in stats.items():
            if proc['state'] == 'Z' or not proc['argv']:
                continue  # zombies are already dead
            try:
                if os.stat(f'/proc/{pid}').st_uid != uid:
                    continue  # never touch other users' processes
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                continue
            name = os.path.basename(proc['argv'][0])
            # scope: the mmseqs workflow of the colabfold MSA runs
            if name != 'mmseqs' and 'blastpgp.sh' not in ' '.join(proc['argv']):
                continue
            if not Clean._in_msa_tree(pid, stats):
                continue  # an unrelated mmseqs run: leave it alone
            if Clean._has_live_driver(pid, stats):
                # running MSA job: keep
                if verbose:
                    print(f"Keep running MSA process: pid={pid} "
                          f"{' '.join(proc['argv'][:2])}")
                kept.append(pid)
                continue
            targets.append(pid)
        if not targets:
            if verbose:
                print('Clean: no unused mmseqs process.')
            return []
        if verbose:
            for pid in targets:
                print(f"Kill unused mmseqs process: pid={pid} "
                      f"{' '.join(stats[pid]['argv'][:2])}")

        # stop the blastpgp.sh wrappers first, so that they cannot
        # spawn the next mmseqs workflow step
        targets.sort(key=lambda pid: 0
                     if 'blastpgp.sh' in ' '.join(stats[pid]['argv']) else 1)

        def _send(pids: list, sig: int) -> list:
            sent = []
            for pid in pids:
                try:
                    os.kill(pid, sig)
                    sent.append(pid)
                except (ProcessLookupError, PermissionError):
                    pass
            return sent

        killed = _send(targets, signal.SIGTERM)
        # wait the grace period, then force kill the survivors
        deadline = time.monotonic() + max(grace, 0.0)
        while time.monotonic() < deadline:
            current = Clean.proc_stats()
            if not any(pid in current and current[pid]['state'] != 'Z'
                       for pid in killed):
                break
            time.sleep(0.2)
        _send(killed, signal.SIGKILL)
        print(f"Clean: killed {len(killed)} unused mmseqs process(es) "
              f"(kept {len(kept)} running MSA process(es))")
        return killed

    @staticmethod
    def clean_pred():
        ''''''
        pass

