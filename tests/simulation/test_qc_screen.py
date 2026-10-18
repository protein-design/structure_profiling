'''
unit tests for qc_screen: backbone_qc_screen() / structure_qc_screen()
- clean model: tests/data/PDB/6BSX_0_A_0.pdb (chain A)
- strained model: tests/data/PDB/2OX8_A.pdb (chain A, vdw clash score 23.19)
- pytest -s tests/simulation/test_qc_screen.py
'''
import os

from ..helper import *
from src.structure_profiling import Pdb
from src.structure_profiling.geometry.geometry import Geometry
from src.structure_profiling.qc.qc_screen import (
    BACKBONE_GATES, STRUCTURE_GATES, backbone_qc_screen, structure_qc_screen)

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')


@ddt
class TestQcScreen(TestCase):

    def test_backbone_screen_clean_pass(self):
        res = backbone_qc_screen(PDB_6BSX, 'A', keep_details=True)
        self.assertTrue(res['pass'])
        self.assertEqual(res['reason'], 'pass')
        self.assertEqual(res['num_missing'], 0)
        self.assertEqual(res['num_breaks'], 0)
        self.assertEqual(res['num_fragments'], 1)
        self.assertEqual(res['num_strained'], 0)
        self.assertAlmostEqual(res['ca_max_dist'], 3.88, places=2)
        self.assertAlmostEqual(res['pep_max_dist'], 1.36, places=2)
        self.assertEqual(res['num_clashes'], 0)
        self.assertEqual(res['num_servere_clashes'], 0)
        self.assertEqual(res['clash_score'], 0.0)

    def test_backbone_rmsz_values(self):
        '''bond_rmsz / angle_rmsz of the screen must match chain_rmsz()'''
        res = backbone_qc_screen(PDB_6BSX, 'A')
        rmsz = Geometry(PDB_6BSX).chain_rmsz(Pdb(PDB_6BSX).get_chain('A'))
        self.assertAlmostEqual(res['bond_rmsz'], rmsz['bond_rmsz'])
        self.assertAlmostEqual(res['angle_rmsz'], rmsz['angle_rmsz'])
        self.assertAlmostEqual(res['bond_rmsz'], 0.6415, places=3)
        self.assertAlmostEqual(res['angle_rmsz'], 0.7344, places=3)
        self.assertEqual(res['n_z_outliers'], 0)
        self.assertEqual(res['z_outlier_frac'], 0.0)
        self.assertEqual(res['peptide_outlier_frac'], 0.0)
        self.assertEqual(res['mainchain_outlier_frac'], 0.0)

    def test_backbone_screen_details(self):
        res = backbone_qc_screen(PDB_6BSX, 'A', keep_details=True)
        self.assertEqual(set(res['details'].keys()),
            {'bond_angle_z', 'planarity', 'chain_breaks', 'chain_strained'})

    def test_structure_screen_clean_pass(self):
        res = structure_qc_screen(PDB_6BSX, 'A', keep_details=True)
        self.assertTrue(res['pass'])
        self.assertEqual(res['reason'], 'pass')
        self.assertEqual(res['cb_outlier_frac'], 0.0)
        self.assertEqual(res['n_chirality_outliers'], 0)
        self.assertAlmostEqual(res['rama_outlier_frac'], 0.006, places=3)
        self.assertAlmostEqual(res['vdw_clash_score'], 3.79, places=2)
        self.assertEqual(set(res['details'].keys()),
            {'cbeta_deviation', 'ca_chirality', 'vdw_clashes', 'rama_outliers',
             'bond_angle_z', 'planarity', 'chain_breaks', 'chain_strained'})

    def test_structure_screen_flags_vdw_clashes(self):
        # 2OX8 side chains overlap: MolProbity score 23.19 > cutoff 20
        res = structure_qc_screen(PDB_2OX8, 'A')
        self.assertFalse(res['pass'])
        self.assertEqual(res['reason'], 'vdw_clashes')
        self.assertAlmostEqual(res['vdw_clash_score'], 23.19, places=2)

    def test_backbone_screen_strained_model_passes(self):
        res = backbone_qc_screen(PDB_2OX8, 'A')
        self.assertTrue(res['pass'])
        self.assertAlmostEqual(res['bond_rmsz'], 0.383, places=3)
        self.assertAlmostEqual(res['angle_rmsz'], 0.982, places=3)
        self.assertEqual(res['num_clashes'], 2)

    def test_no_chain(self):
        res = backbone_qc_screen(PDB_6BSX, 'Z')
        self.assertFalse(res['pass'])
        self.assertEqual(res['reason'], 'no_chain')

    @data(['ca_max_dist', 'ca_distance'],
          ['pep_max_dist', 'peptide_distance'],
          ['bond_rmsz', 'bond_rmsz'],
          ['angle_rmsz', 'angle_rmsz'])
    @unpack
    def test_gate_override_fails(self, gate, expect_reason):
        # zero tolerance on one gate -> exactly that gate must fail
        res = backbone_qc_screen(PDB_6BSX, 'A', gates={gate: 0.0})
        self.assertFalse(res['pass'])
        self.assertEqual(res['reason'], expect_reason)

    def test_no_details_by_default(self):
        res = backbone_qc_screen(PDB_6BSX, 'A', keep_details=False)
        self.assertIsNone(res['details'])

    def test_default_gates(self):
        self.assertEqual(BACKBONE_GATES['num_missing'], 0)
        self.assertEqual(BACKBONE_GATES['num_breaks'], 0)
        self.assertEqual(BACKBONE_GATES['num_fragments'], 1)
        self.assertEqual(BACKBONE_GATES['bond_rmsz'], 2.0)
        self.assertEqual(BACKBONE_GATES['angle_rmsz'], 1.5)
        self.assertEqual(BACKBONE_GATES['ca_max_dist'], 5.0)
        self.assertEqual(BACKBONE_GATES['pep_max_dist'], 5.0)
        # structure gates extend the backbone ones
        for key, value in BACKBONE_GATES.items():
            self.assertEqual(STRUCTURE_GATES[key], value)
        self.assertEqual(STRUCTURE_GATES['cb_outlier_frac'], 0.05)
        self.assertEqual(STRUCTURE_GATES['rama_outlier_frac'], 0.05)
        self.assertEqual(STRUCTURE_GATES['n_chirality_outliers'], 0)
        self.assertEqual(STRUCTURE_GATES['vdw_clash_score'], 20.0)
