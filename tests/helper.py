from copy import deepcopy
from datetime import datetime
import json
import os
import pandas as pd
import shutil
import sys
from pprint import pprint
from ddt import ddt, data, unpack
from unittest import TestCase, mock, skip

# store example data for unit and local testing
TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, 'data')
RESULTS_DIR = os.path.join(TESTS_DIR, 'results')
# legacy alias still referenced by several test modules
DIR_DATA = DATA_DIR

#
profiling_dir = os.path.dirname(TESTS_DIR)
src_dir = os.path.join(profiling_dir, 'src')

from dotenv import load_dotenv
env_file = os.path.join(profiling_dir, '.env')
load_dotenv(env_file)
