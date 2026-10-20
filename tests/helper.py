from copy import deepcopy
from datetime import datetime
import json
import os
import numpy as np
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


# ---------------------------------------------------------------------------
# PDB editing: synthesize chain breaks / strain / missing atoms on the fly
# (all PDB test data is clean, so defect models are written to temp files)
# ---------------------------------------------------------------------------

def pdb_coord(line):
    '''(x, y, z) of an ATOM/HETATM record (PDB fixed columns 31-54).'''
    return np.array([float(line[30:38]), float(line[38:46]),
                     float(line[46:54])])


def set_pdb_coord(line, coord):
    '''ATOM/HETATM record with the coordinate columns 31-54 replaced.'''
    return line[:30] + f'{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}' \
        + line[54:]


def _is_atom(line):
    return line.startswith(('ATOM', 'HETATM'))


def _resseq(line):
    return int(line[22:26])


def _atom_name(line):
    return line[12:16].strip()


def write_pdb_lines(lines, pdb_file):
    '''write PDB text lines back to a file (returns pdb_file).'''
    with open(pdb_file, 'w') as f:
        f.writelines(lines)
    return pdb_file


def translate_residues(lines, start_resseq, offset):
    '''rigid-body translate every ATOM/HETATM record with
    resseq >= start_resseq by offset (A).

    A large offset (e.g. (20, 0, 0)) splits the chain into two rigid
    fragments: exactly one break at the (start_resseq-1, start_resseq)
    peptide, all remaining pairs intact.
    '''
    return [
        set_pdb_coord(l, pdb_coord(l) + np.asarray(offset, dtype=float))
        if _is_atom(l) and _resseq(l) >= start_resseq else l
        for l in lines
    ]


def delete_residue(lines, resseq):
    '''drop every ATOM/HETATM record of one residue.

    Creates a residue-numbering gap (resnum_gap = 2) on top of the
    broken peptide between resseq-1 and resseq+1.
    '''
    return [l for l in lines if not (_is_atom(l) and _resseq(l) == resseq)]


def strip_atom(lines, resseq, atom_name):
    '''drop one named atom record of one residue (e.g. 'O').

    Leaves the residue with missing backbone atoms (num_missing).
    '''
    return [l for l in lines
            if not (_is_atom(l) and _resseq(l) == resseq
                    and _atom_name(l) == atom_name)]


def stretch_peptide(lines, resseq, extra=0.5):
    '''move N of resseq `extra` A along the C(resseq-1) -> N(resseq)
    direction.

    Stretches exactly that peptide bond into the strained window
    (Breaks.cn_strain < d_C_N <= Breaks.cn_break, e.g. 1.33 + 0.5 =
    1.83 A) without breaking it; CA-CA and all other pairs unchanged.
    '''
    c_prev = n_next = None
    for l in lines:
        if _is_atom(l):
            if _resseq(l) == resseq - 1 and _atom_name(l) == 'C':
                c_prev = pdb_coord(l)
            if _resseq(l) == resseq and _atom_name(l) == 'N':
                n_next = pdb_coord(l)
    if c_prev is None or n_next is None:
        raise ValueError(f'peptide {resseq - 1}-{resseq} not found')
    direction = (n_next - c_prev) / np.linalg.norm(n_next - c_prev)
    return [
        set_pdb_coord(l, pdb_coord(l) + extra * direction)
        if _is_atom(l) and _resseq(l) == resseq and _atom_name(l) == 'N'
        else l
        for l in lines
    ]
