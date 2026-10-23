'''
unit tests for clash analysis (structure/clashes.py, Clashes):
steric clashes, severe clashes, clash score, MolProbity vdW overlaps
and the run() / check() summary + gates
- clean model: tests/data/PDB/6BSX_0_A_0.pdb
    (chain A: 170 residues, 1319 heavy atoms; no steric clashes,
     5 vdW overlaps, clashscore 3.79 -> pass)
- strained model: tests/data/PDB/2OX8_A.pdb
    (chain A: 1078 heavy atoms; 2 steric clashes, none severe (< 2.0 A),
     25 vdW overlaps, clashscore 23.19 -> vdw_clashes gate)
- pytest -s tests/structure/test_clashes.py
'''
import os

import numpy as np
import pandas as pd

from ..helper import *
from src.structure_profiling import Pdb, Clashes
from src.structure_profiling.geometry.geometry_distance import GeometryDistance

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')


class TestClashes(TestCase):

    def setUp(self):
        self.chain = Pdb(PDB_6BSX).get_chain('A')
        self.clashes = Clashes(self.chain)

    def test_standard_residues(self):
        self.assertEqual(len(self.clashes.standard_residues()), 170)

    def test_heavy_atoms(self):
        self.assertEqual(len(self.clashes.heavy_atoms()), 1319)

    def test_steric_clashes(self):
        # clean model: no clashes at all
        steric = self.clashes.steric_clashes()
        self.assertEqual(len(steric), 0)
        self.assertEqual(list(steric.columns),
                         ['atom_pos', 'atom', 'other_pos', 'other_atom',
                          'distance'])
        # strained model: two contacts below 2.4 A, none below 2.0 A
        strained = Clashes(Pdb(PDB_2OX8).get_chain('A'))
        steric = strained.steric_clashes()
        self.assertEqual(len(steric), 2)
        self.assertGreaterEqual(steric['distance'].min(), 2.0)
        # a wider cutoff can only add pairs
        self.assertGreaterEqual(len(strained.steric_clashes(cutoff=3.0)), 2)

    def test_severe_clashes(self):
        strained = Clashes(Pdb(PDB_2OX8).get_chain('A'))
        steric = strained.steric_clashes()
        # all contacts >= 2.0 A: not severe
        self.assertEqual(len(strained.severe_clashes(steric)), 0)
        # below the clash cutoff everything counts severe
        self.assertEqual(len(strained.severe_clashes(steric, distance=2.4)), 2)
        # empty table passes through unchanged
        self.assertTrue(strained.severe_clashes(steric.iloc[0:0]).empty)

    def test_clash_score(self):
        strained = Clashes(Pdb(PDB_2OX8).get_chain('A'))
        steric = strained.steric_clashes()
        score = strained.clash_score(steric)
        self.assertEqual(score, 1.86)
        self.assertEqual(score, round(1000.0 * 2 / 1078, 2))
        # explicit n_atoms is honoured, no atoms -> 0
        self.assertEqual(strained.clash_score(steric, n_atoms=1000), 2.0)
        self.assertEqual(self.clashes.clash_score(steric, n_atoms=0), 0.0)

    def test_vdw_clashes(self):
        vdw = self.clashes.vdw_clashes(include_hydrogen=False)
        self.assertEqual(len(vdw), 5)
        self.assertEqual(list(vdw.columns),
                         ['atom_pos', 'atom', 'element', 'other_pos',
                          'other_atom', 'other_element', 'distance',
                          'overlap'])
        # MolProbity: overlap = (r_i + r_j) - d_ij >= 0.40 A
        self.assertGreaterEqual(vdw['overlap'].min(), Clashes.VDW_OVERLAP_CUT)
        # recompute the first overlap from the recorded distance
        row = vdw.iloc[0]
        radii = Clashes.VDW_RADII
        expect = radii[row['element']] + radii[row['other_element']] \
            - row['distance']
        self.assertAlmostEqual(float(row['overlap']), expect, places=6)
        # strained model: many side-chain overlaps
        strained = Clashes(Pdb(PDB_2OX8).get_chain('A'))
        self.assertEqual(len(strained.vdw_clashes(include_hydrogen=False)), 25)

    def test_vdw_clash_score(self):
        vdw = self.clashes.vdw_clashes(include_hydrogen=False)
        self.assertEqual(self.clashes.vdw_clash_score(vdw),
                         {'clash_score': 3.79, 'n_clashes': 5,
                          'max_overlap': 0.609})
        # no clashes / no atoms -> zeros
        zeros = {'clash_score': 0.0, 'n_clashes': 0, 'max_overlap': 0.0}
        self.assertEqual(self.clashes.vdw_clash_score(vdw.iloc[0:0]), zeros)
        self.assertEqual(self.clashes.vdw_clash_score(vdw, n_atoms=0), zeros)

    def test_run_summary(self):
        info = self.clashes.run(keep_details=False)
        self.assertEqual(set(info), {'n_heavy_atoms', 'num_clashes',
            'num_severe_clashes', 'clash_score', 'n_vdw_clashes',
            'vdw_clash_score', 'max_overlap'})
        self.assertEqual(info['n_heavy_atoms'], 1319)
        self.assertEqual(info['num_clashes'], 0)
        self.assertEqual(info['num_severe_clashes'], 0)
        self.assertEqual(info['clash_score'], 0.0)
        self.assertEqual(info['n_vdw_clashes'], 5)
        self.assertEqual(info['vdw_clash_score'], 3.79)
        # details hold the raw tables
        detailed = self.clashes.run()
        self.assertIn('steric_clashes', detailed['details'])
        self.assertIn('vdw_clashes', detailed['details'])

    def test_check(self):
        # clean model passes all gates
        self.assertEqual(self.clashes.check(), (True, 'pass'))
        # strained model: MolProbity clashscore above the poor band
        strained = Clashes(Pdb(PDB_2OX8).get_chain('A'))
        self.assertEqual(strained.check(), (False, 'vdw_clashes'))
        # synthetic gates: severe clash, steric score, null summary
        self.assertEqual(strained.check({'num_severe_clashes': 1}),
                         (False, 'severe_clashes'))
        self.assertEqual(strained.check({'clash_score': 31.0}),
                         (False, 'clashes'))
        self.assertEqual(strained.check({}), (False, 'null_clashes'))

    def test_geometry_delegates_agree(self):
        # GeometryDistance clash APIs must stay consistent with Clashes
        gd = GeometryDistance(self.chain)
        self.assertTrue(gd.steric_clashes().equals(
            self.clashes.steric_clashes()))
        self.assertTrue(gd.vdw_clashes().equals(self.clashes.vdw_clashes()))
        vdw = self.clashes.vdw_clashes()
        self.assertEqual(gd.vdw_clash_score(vdw),
                         self.clashes.vdw_clash_score(vdw))

    def test_cutoffs_shared_with_geometry(self):
        from src.structure_profiling.geometry.geometry import Geometry
        self.assertEqual(Geometry.clash_score_cutoff,
                         Clashes.clash_score_cutoff)
        self.assertEqual(Geometry.severe_clash_distance,
                         Clashes.severe_clash_distance)
        self.assertEqual(Geometry.vdw_clash_cutoff, Clashes.vdw_clash_cutoff)
