#!/usr/bin/python
'''
requirements:
    - none: read-only diagnostics of the runtime environment,
      no table and no pdb_group involved
example:
    - sp test_env
    - python -m src.pipelines.test_env
function:
    - report python version and imports/versions of every dependency
      (required: numpy, pandas, biopython, mysql-connector-python,
      python-dotenv, openmm, pdbfixer; optional: torch, scipy, matplotlib, mdtraj)
    - check the .env resolved by src/helper.py: path, loaded keys (secrets
      masked), required keys, output dirs exist and are writable
    - check accelerators: OpenMM platforms (CUDA availability), nvidia-smi
      GPUs, PyTorch CUDA when torch is installed
    - check external binaries (mkdssp, needed by sp chain_dssp)
    - probe MySQL with the .env credentials
exit:
    - status 1 when a required check fails (MISS), 0 otherwise (warnings OK)
'''
import os
import shutil
import subprocess
import sys
from importlib import import_module, metadata
from pathlib import Path

from dotenv import dotenv_values
from openmm import Platform

from src.helper import *

# modules every pipeline depends on
REQUIRED_MODULES = ('numpy', 'pandas', 'Bio', 'mysql.connector', 'dotenv', 'openmm', 'pdbfixer')
# nice to have, only reported when installed
OPTIONAL_MODULES = ('torch', 'scipy', 'matplotlib', 'mdtraj')
# pip distribution names differing from the import name
DIST_NAMES = {
    'Bio': 'biopython',
    'mysql.connector': 'mysql-connector-python',
    'dotenv': 'python-dotenv',
}
# keys the database layer (structure_profiling/database/db.py) and helper need
REQUIRED_ENV_KEYS = ('MYSQL_HOST', 'MYSQL_USER', 'MYSQL_PASSWORD', 'MYSQL_DB', 'output_dir')
# substrings marking a secret .env value that must never be printed
SECRET_MARKS = ('PASSWORD', 'PASSWD', 'SECRET', 'TOKEN')
# dirs configured in helper params, checked for existence and write access
CHECKED_DIRS = ('output_dir', 'output_pdb_dir', 'tmp_dir', 'rawdata_dir')
# mkdssp 4.x (libcifpp based) needs the Chemical Component Dictionary at
# runtime, otherwise every sp chain_dssp run aborts with a "Configuration
# error" about a missing components.cif. First hit wins; note libcifpp
# also honours the LIBCIFPP_DATA_DIR environment variable on its own.
CCD_NAMES = ('components.cif', 'components.cif.gz')
CCD_DIRS = (
    os.getenv('LIBCIFPP_DATA_DIR'),
    Path.home() / '.local/share/libcifpp',
    Path('/var/cache/libcifpp'),
    Path('/usr/share/libcifpp'),
    Path('/usr/local/share/libcifpp'),
)

OK, WARN, MISS = 'OK', 'WARN', 'MISS'


def module_version(name:str):
    '''best-effort version string of an importable module'''
    mod = import_module(name)
    version = getattr(mod, '__version__', None)
    if version is None:
        try:
            version = metadata.version(DIST_NAMES.get(name, name))
        except Exception:
            version = 'unknown'
    return str(version)


def check_imports():
    '''import python and every dependency, one report row per module'''
    rows = [{
        'item': 'python',
        'status': OK if sys.version_info >= (3, 11) else MISS,
        'detail': f"{sys.version.split()[0]} ({sys.executable})",
    }]
    for name in REQUIRED_MODULES + OPTIONAL_MODULES:
        required = name in REQUIRED_MODULES
        try:
            rows.append({'item': name, 'status': OK, 'detail': module_version(name)})
        except Exception as e:
            rows.append({
                'item': name,
                'status': MISS if required else WARN,
                'detail': f'{type(e).__name__}: {e}',
            })
    return rows


def check_env_file():
    '''.env resolved by src/helper.py: presence, keys (secrets masked), dirs'''
    rows = []
    env_path = Path(env_file)
    rows.append({
        'item': '.env path',
        'status': OK if env_path.is_file() else MISS,
        'detail': str(env_path),
    })
    if not env_path.is_file():
        return rows

    # keys present in the file (secrets never printed)
    values = dotenv_values(env_path)
    loaded = sorted(k for k, v in values.items() if k and v is not None)
    rows.append({
        'item': 'keys loaded',
        'status': OK if loaded else MISS,
        'detail': ', '.join(loaded) or 'none',
    })
    missing = [k for k in REQUIRED_ENV_KEYS if not os.getenv(k)]
    rows.append({
        'item': 'required keys',
        'status': MISS if missing else OK,
        'detail': ', '.join(missing) if missing else f"all {len(REQUIRED_ENV_KEYS)} present",
    })
    for k in loaded:
        if any(s in k.upper() for s in SECRET_MARKS) or k in CHECKED_DIRS:
            continue
        v = os.getenv(k)
        if v:
            rows.append({'item': k, 'status': OK, 'detail': v})

    # output dirs from helper params
    for name in CHECKED_DIRS:
        path = params.get(name)
        if path is None:
            continue
        path = Path(path)
        if not path.is_dir():
            rows.append({'item': name, 'status': MISS, 'detail': f'{path} (not a directory)'})
        elif not os.access(path, os.W_OK):
            rows.append({'item': name, 'status': MISS, 'detail': f'{path} (not writable)'})
        else:
            rows.append({'item': name, 'status': OK, 'detail': str(path)})
    return rows


def check_cuda():
    '''OpenMM platforms, driver GPUs and PyTorch CUDA (torch optional)'''
    rows = []
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        speed = ''
        try:
            speed = f' speed={p.getSpeedEstimate():.0f}'
        except Exception:
            pass
        rows.append({'item': f'openmm platform {i}', 'status': OK, 'detail': f'{p.getName()}{speed}'})

    try:
        cuda = Platform.getPlatformByName('CUDA')
        detail = 'platform registered'
        try:
            detail += f", default precision={cuda.getPropertyDefaultValue('CudaPrecision')}"
        except Exception:
            pass
        rows.append({'item': 'openmm CUDA', 'status': OK, 'detail': detail})
    except Exception as e:
        rows.append({'item': 'openmm CUDA', 'status': WARN, 'detail': f'not available: {e}'})

    # GPUs reported by the driver
    smi = shutil.which('nvidia-smi')
    if smi is None:
        rows.append({'item': 'nvidia-smi', 'status': WARN, 'detail': 'binary not on PATH'})
    else:
        try:
            res = subprocess.run(
                [smi, '--query-gpu=name,memory.total,driver_version', '--format=csv,noheader'],
                capture_output=True, text=True, timeout=15,
            )
            gpus = [line for line in res.stdout.strip().splitlines() if line]
            rows.append({
                'item': 'nvidia-smi',
                'status': OK if res.returncode == 0 and gpus else WARN,
                'detail': ' | '.join(gpus) or res.stderr.strip(),
            })
        except Exception as e:
            rows.append({'item': 'nvidia-smi', 'status': WARN, 'detail': f'{type(e).__name__}: {e}'})
    visible = os.getenv('CUDA_VISIBLE_DEVICES')
    rows.append({
        'item': 'CUDA_VISIBLE_DEVICES',
        'status': OK,
        'detail': visible if visible else 'not set (all GPUs visible)',
    })

    # PyTorch CUDA when torch happens to be installed
    try:
        import torch
        detail = f'torch {torch.__version__}, cuda {torch.version.cuda}, available={torch.cuda.is_available()}'
        if torch.cuda.is_available():
            names = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
            detail += f', devices={names}'
        rows.append({'item': 'pytorch cuda', 'status': OK, 'detail': detail})
    except Exception as e:
        rows.append({'item': 'pytorch cuda', 'status': WARN, 'detail': f'torch not installed ({type(e).__name__})'})
    return rows


def find_ccd():
    '''path of the Chemical Component Dictionary file mkdssp 4.x needs'''
    for data_dir in CCD_DIRS:
        if not data_dir:
            continue
        for name in CCD_NAMES:
            path = Path(data_dir) / name
            if path.is_file():
                return path
    return None


def check_binaries():
    '''external programs some pipelines shell out to'''
    rows = []
    mkdssp = shutil.which('mkdssp')
    rows.append({
        'item': 'mkdssp',
        'status': OK if mkdssp else WARN,
        'detail': mkdssp or 'not on PATH (sp chain_dssp will fail)',
    })
    if mkdssp:
        ccd = find_ccd()
        if ccd is not None:
            rows.append({'item': 'mkdssp CCD', 'status': OK, 'detail': str(ccd)})
        else:
            fix = (
                'curl -o ~/.local/share/libcifpp/components.cif.gz '
                'https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz'
            )
            rows.append({
                'item': 'mkdssp CCD',
                'status': WARN,
                'detail': f'components.cif missing (mkdssp 4.x aborts without it); fix: {fix}',
            })
    return rows


def check_mysql():
    '''probe the database with the .env credentials (60 s connect timeout)'''
    try:
        rows = QueryComplex(False).list_data('SELECT VERSION() AS server_version, DATABASE() AS db')
        info = rows[0] if rows else {}
        return [{
            'item': 'mysql',
            'status': OK,
            'detail': f"server {info.get('server_version')}, db {info.get('db')}",
        }]
    except Exception as e:
        return [{'item': 'mysql', 'status': WARN, 'detail': f'{type(e).__name__}: {e}'}]


def report(title, rows):
    '''print one section of the environment report'''
    print(f'\n### {title}')
    for row in rows:
        print(f"  [{row['status']:^4}] {row['item']:22s} {row['detail']}")
    return rows


def run(params, meta):
    '''execute all checks, print the report and tally statuses into meta'''
    sections = (
        ('python & dependencies', check_imports()),
        ('.env file', check_env_file()),
        ('cuda & accelerators', check_cuda()),
        ('external binaries', check_binaries()),
        ('mysql connectivity', check_mysql()),
    )
    for title, rows in sections:
        report(title, rows)
        for row in rows:
            meta[row['status']] += 1
    return sections


if __name__ == '__main__':
    run(params, meta)
    print(f"\nSUMMARY: {meta[OK]} ok, {meta[WARN]} warn, {meta[MISS]} miss")
    footer(meta, params)
    if meta[MISS]:
        sys.exit(1)