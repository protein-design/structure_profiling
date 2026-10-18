"""
IO intensive writing into Linux temporary directory
"""
import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import hashlib
import subprocess
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class WriteTmp:
    def __init__(self, temporary_dir:str, destination_dir:str, name:str=None):
        """
        destination_dir: final directory
        temporary_dir: temporary directory
        """
        self.temp_base = Path(temporary_dir)
        # NB: attribute was misspelled 'destinaton_dir'
        self.destination_dir = Path(destination_dir)
        self.destination_dir.mkdir(parents=True, exist_ok=True)
        self.name = name if name else 'io'
        # NB: initialize so that temp_path raises a clean RuntimeError
        # (AttributeError) before __enter__ was called
        self.temp_dir_path = None

    def __enter__(self):
        # create temp_path
        user_id = os.getuid()
        pid = os.getpid()
        prefix = f"{self.name}_{user_id}_{pid}_"
        self.temp_dir_path = Path(tempfile.mkdtemp(
            prefix=prefix,
            dir=self.temp_base
        ))
        logger.info(f"create temporary directory: {self.temp_dir_path}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # move to destination
        self.copy_temp()

        # delete temp_path
        if self.temp_dir_path and self.temp_dir_path.exists():
            try:
                shutil.rmtree(self.temp_dir_path)
                logger.info(f"clean temporary directory: {self.temp_dir_path}")
            except Exception as e:
                logger.warning(f"fail: err={e}")

    @property
    def temp_path(self) -> Path:
        if not self.temp_dir_path:
            raise RuntimeError(f"{self.temp_dir_path} is not existing.")
        return self.temp_dir_path
    
    def copy_temp(self):
        if not self.temp_dir_path:
            return False
        # move files rather than copy
        cmd = [
            'rsync', '-av', '--progress', '--remove-source-files',
            str(self.temp_dir_path) + '/',
            str(self.destination_dir) + '/',
        ]
        print('Move results from tempory to destination:', cmd)
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"rsync failed: {result.stderr}")
            return False
        logger.info(f"rsync is done: {result.stdout}")
        return True