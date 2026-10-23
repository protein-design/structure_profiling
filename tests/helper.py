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


def move_atom(lines, resseq, atom_name, coord):
    '''set the coordinate of one named atom of one residue to `coord`.

    A radial move (e.g. N pulled along its own CA->N direction out to
    4 A) keeps the intra-residue N-CA-C angle but deforms everything
    around the atom: the bond angles measured at / across it collapse
    (|Z| >> 5) and its peptide bond is strained or broken.
    '''
    coord = np.asarray(coord, dtype=float)
    return [
        set_pdb_coord(l, coord)
        if _is_atom(l) and _resseq(l) == resseq
        and _atom_name(l) == atom_name else l
        for l in lines
    ]


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


def rotate_residues(lines, start_resseq, angle=180.0):
    '''rigid-body rotate each residue with resseq >= start_resseq in
    place around its own CA atom, by `angle` degrees (rotation about x
    then y, determinant +1).

    Per-residue rotation: intra-residue bond angles, bond lengths and
    CA chirality are preserved (rigid proper rotation) while everything
    crossing residues downstream of peptide start_resseq-1 ->
    start_resseq is scrambled - the peptide-bridge bond angles
    (CA-C-N, C-N-CA) leave the Engh & Huber window, omega leaves the
    cis/trans windows (Torsions.OMEGA_DEV_CUT = 30 deg about 0 / 180)
    and phi / psi become Ramachandran garbage. Use it to synthesize
    torsion / planarity / covalent-geometry outliers without breaking
    residues apart.
    '''
    a = np.radians(angle)
    rx = np.array([[1, 0, 0],
                   [0, np.cos(a), -np.sin(a)],
                   [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(a), 0, np.sin(a)],
                   [0, 1, 0],
                   [-np.sin(a), 0, np.cos(a)]])
    rot = rx @ ry
    pivots = {}
    for l in lines:
        if _is_atom(l) and _resseq(l) >= start_resseq \
                and _atom_name(l) == 'CA':
            pivots[_resseq(l)] = pdb_coord(l)
    if not pivots:
        raise ValueError(f'no residues with resseq >= {start_resseq}')
    return [
        set_pdb_coord(l, rot @ (pdb_coord(l) - pivots[_resseq(l)])
                      + pivots[_resseq(l)])
        if _is_atom(l) and _resseq(l) >= start_resseq else l
        for l in lines
    ]

# ---------------------------------------------------------------------------
# synthetic PDB with a HETATM ligand (the shipped test PDBs have no HETATM)
# ---------------------------------------------------------------------------

def atom_record(rec_name, serial, atom_name, resname, chain_no, resseq,
                coord, element):
    '''one ATOM / HETATM record with fixed-format columns (PDB v3).'''
    return (f'{rec_name:<6}{serial:>5} {atom_name:>4}'
            f'{" ":1}{resname:>3} {chain_no:1}{resseq:>4}{" ":1}   '
            f'{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}'
            f'{1.00:6.2f}{0.00:6.2f}          {element:>2}  \n')


def synthesize_residue(serial_start, resseq, coord, z=0.0, chain_no='A',
                       resname='ALA'):
    '''five heavy atoms (N, CA, C, O, CB) of one residue placed around
    coord = (x, y) at height z: N (x-1.3, 0.8), CA (x, 0), C (x+1.3,
    -0.6), O (x+1.3, -1.8), CB (x, 1.5). Returns (lines, next_serial).
    '''
    x, y = float(coord[0]), float(coord[1])
    atoms = [('N', (x - 1.3, y + 0.8), 'N'), ('CA', (x, y), 'C'),
             ('C', (x + 1.3, y - 0.6), 'C'), ('O', (x + 1.3, y - 1.8), 'O'),
             ('CB', (x, y + 1.5), 'C')]
    lines, serial = [], serial_start
    for name, c, element in atoms:
        lines.append(atom_record('ATOM', serial, name, resname, chain_no,
                                 resseq, (c[0], c[1], z), element))
        serial += 1
    return lines, serial


def synthesize_site_pdb(pdb_file):
    '''small synthetic chain with a HETATM ligand + water, exact coords:

    chain A holds two sequence fragments (residues 1-3 and 101-103 at
    z = 0 and z = 8, 3.8 A CA spacing along x) plus ligand LIG 300
    (C1 (3.8, 2.0, 4.0), N1 (3.8, 3.4, 4.0)) and water HOH 501
    ((7.6, 2.5, 4.0)). At the default 5.0 A cutoff the ligand touches
    residues 2 and 102 (min CB-C1 sqrt(16.25) ~ 4.031 A, 4 contact
    atoms each) and 3 / 103 (N, 4.868 A): neighbors close in space but
    scattered over the sequence (n_fragments = 2). Returns pdb_file.
    '''
    lines, serial = [], 1
    for resseq, z in ((1, 0.0), (2, 0.0), (3, 0.0),
                      (101, 8.0), (102, 8.0), (103, 8.0)):
        res_lines, serial = synthesize_residue(
            serial, resseq, (3.8 * ((resseq - 1) % 100), 0.0), z=z)
        lines += res_lines
    # ligand LIG 300 (hetflag H_LIG) and water HOH 501 (hetflag W)
    lines.append(atom_record('HETATM', serial, 'C1', 'LIG', 'A', 300,
                             (3.8, 2.0, 4.0), 'C'))
    lines.append(atom_record('HETATM', serial + 1, 'N1', 'LIG', 'A', 300,
                             (3.8, 3.4, 4.0), 'N'))
    lines.append(atom_record('HETATM', serial + 2, 'O', 'HOH', 'A', 501,
                             (7.6, 2.5, 4.0), 'O'))
    lines.append('END\n')
    return write_pdb_lines(lines, pdb_file)

