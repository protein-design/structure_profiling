'''
unit tests for protein geometry QC:
GeometryDistance / GeometryAngle / RamaGrid / Geometry
- clean model: tests/data/PDB/6BSX_0_A_0.pdb
    (chain A: 170 residues, numbering starts at 2; no breaks, no clashes)
- strained model: tests/data/PDB/2OX8_A.pdb
    (chain A: 129 residues; 2 steric clashes, none severe (< 2.0 A))
- pytest -s tests/simulation/test_geometry.py
'''
import math
import os
from unittest import skipUnless

import numpy as np
import pandas as pd

from ..helper import *
from src.structure_profiling import Pdb
from src.structure_profiling.qc.geometry import Geometry
from src.structure_profiling.qc.geometry_angle import GeometryAngle
from src.structure_profiling.qc.geometry_distance import GeometryDistance
from src.structure_profiling.qc.rama_grid import RamaGrid

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')
RAMA_NPZ = os.path.join(src_dir, 'structure_profiling', 'data', 'Rama8000',
                        'rama8000-general-grids.npz')


class TestGeometryDistance(TestCase):

    def setUp(self):
        self.chain = Pdb(PDB_6BSX).get_chain('A')
        self.gd = GeometryDistance(self.chain)

    def test_standard_residues(self):
        self.assertEqual(len(self.gd.standard_residues()), 170)

    def test_heavy_atoms(self):
        self.assertEqual(len(self.gd.heavy_atoms()), 1319)

    def test_atom_distance(self):
        # consecutive CA atoms (6BSX chain A numbering starts at residue 2)
        self.assertAlmostEqual(self.gd.atom_distance(2, 'CA', 3, 'CA'),
                               3.792, places=3)
        # missing residue / atom -> None
        self.assertIsNone(self.gd.atom_distance(99999, 'CA', 3, 'CA'))
        self.assertIsNone(self.gd.atom_distance(2, 'XX', 3, 'CA'))

    def test_ca_and_peptide_distances(self):
        ca = self.gd.ca_distances()
        pep = self.gd.peptide_distances()
        self.assertAlmostEqual(float(np.mean(ca)), 3.805, places=3)
        self.assertAlmostEqual(float(np.mean(pep)), 1.331, places=3)

    def test_steric_clashes(self):
        # clean model: no clashes at all
        self.assertEqual(len(self.gd.steric_clashes()), 0)
        # strained model: two contacts below 2.4 A, none below 2.0 A
        gd = GeometryDistance(Pdb(PDB_2OX8).get_chain('A'))
        clashes = gd.steric_clashes()
        self.assertEqual(len(clashes), 2)
        self.assertEqual(list(clashes.columns),
            ['atom_pos', 'atom', 'other_pos', 'other_atom', 'distance'])
        self.assertTrue((clashes['distance'] < 2.4).all())
        self.assertGreaterEqual(clashes['distance'].min(), 2.0)

    def test_bond_length_z(self):
        z = self.gd.bond_length_z()
        self.assertEqual(list(z.columns),
            ['res_pos', 'resname', 'atoms', 'observed', 'ideal',
             'sd', 'z', 'outlier'])
        n_ca = z[z['atoms'] == 'N-CA']
        self.assertEqual(len(n_ca), 170)
        # observed N-CA matches the Engh & Huber ideal (1.458 A)
        self.assertAlmostEqual(n_ca['observed'].mean(), 1.458, places=2)
        self.assertAlmostEqual(n_ca['ideal'].iloc[0], 1.458, places=3)
        # a refined structure has no covalent outliers
        self.assertFalse(z['outlier'].any())

    def test_cbeta_deviation(self):
        cb = self.gd.cbeta_deviation()
        self.assertEqual(len(cb), 170)
        self.assertEqual(int(cb['outlier'].sum()), 0)
        self.assertLess(cb['cb_dev'].dropna().mean(), 0.1)

    def test_vdw_clashes_and_score(self):
        clashes = self.gd.vdw_clashes(include_hydrogen=False)
        self.assertEqual(len(clashes), 5)
        self.assertGreaterEqual(clashes['overlap'].min(), 0.40)
        score = self.gd.vdw_clash_score(clashes)
        self.assertAlmostEqual(score['clash_score'], 3.79, places=2)
        self.assertEqual(score['n_clashes'], 5)
        self.assertAlmostEqual(score['max_overlap'], 0.609, places=3)


@ddt
class TestGeometryAngle(TestCase):

    def setUp(self):
        self.chain = Pdb(PDB_6BSX).get_chain('A')
        self.ga = GeometryAngle(self.chain)

    @data(
        [90.0, [1, 0, 0], [0, 0, 0], [0, 1, 0]],
        [60.0, [1, 0, 0], [0, 0, 0], [0.5, math.sqrt(3) / 2, 0]],
    )
    @unpack
    def test_bond_angle(self, expect, a, b, c):
        self.assertAlmostEqual(
            GeometryAngle.bond_angle(np.array(a), np.array(b), np.array(c)),
            expect)

    @data(
        [0.0, [0, 1, 1], [0, 0, 0], [1, 0, 0], [1, 1, 1]],
        [90.0, [0, 1, 1], [0, 0, 0], [1, 0, 0], [1, -1, 1]],
        [180.0, [0, 1, 1], [0, 0, 0], [1, 0, 0], [1, -1, -1]],
    )
    @unpack
    def test_dihedral_angle(self, expect, p0, p1, p2, p3):
        self.assertAlmostEqual(
            GeometryAngle.dihedral_angle(np.array(p0), np.array(p1),
                                         np.array(p2), np.array(p3)),
            expect)

    @data(
        [-60, -45, 'right-handed-helix'],
        [-135, 140, 'beta-sheet'],
        [60, 40, 'left-handed-helix'],
        [-45, 120, 'allowed'],
        [60, 170, 'bad'],
        [None, None, 'null'],
        [float('nan'), 0, 'null'],
    )
    @unpack
    def test_evaluate_ramachandran(self, phi, psi, expect):
        self.assertEqual(GeometryAngle.evaluate_ramachandran(phi, psi), expect)

    @data(
        # general residue: alpha / beta favoured, stretched outlier
        ['ALA', -62, -41, None, 'favoured'],
        ['ALA', -135, 140, None, 'favoured'],
        # inside the general beta favoured box
        ['ALA', -45, 120, None, 'favoured'],
        ['ALA', 150, -150, None, 'outlier'],
        # glycine has its own (wider) boxes
        ['GLY', 60, 40, None, 'allowed'],
        # proline: phi restricted near -60
        ['PRO', -60, 140, None, 'favoured'],
        ['PRO', -100, 0, None, 'allowed'],
        ['PRO', 0, 0, None, 'outlier'],
        # pre-proline extension
        ['ALA', -100, 150, 'PRO', 'favoured'],
        # missing torsions
        ['ALA', None, None, None, 'NA'],
    )
    @unpack
    def test_rama_class(self, resname, phi, psi, prev, expect):
        self.assertEqual(GeometryAngle.rama_class(phi, psi, resname, prev),
                         expect)

    def test_bond_angles(self):
        ba = self.ga.bond_angles()
        self.assertEqual(len(ba), 170)
        self.assertAlmostEqual(ba['angle_N_CA_C'].mean(), 111.056, places=2)
        self.assertEqual(int(ba['angle_CA_C_N'].notna().sum()), 169)

    def test_torsion_angles(self):
        ta = self.ga.torsion_angles()
        self.assertEqual(len(ta), 170)
        self.assertTrue(math.isnan(ta.iloc[0]['phi']))  # N-terminal phi
        self.assertAlmostEqual(ta.iloc[0]['omega'], -178.99, places=2)
        self.assertEqual(ta['conformation'].value_counts().to_dict(),
                         {'trans': 169, 'null': 1})

    def test_ca_chirality(self):
        ch = self.ga.ca_chirality()
        self.assertEqual(ch['handedness'].value_counts().to_dict(),
                         {'L': 160, 'NA': 10})   # NA = 10 glycines

    def test_peptide_planarity(self):
        pl = self.ga.peptide_planarity()
        self.assertEqual(len(pl), 169)
        self.assertEqual(int(pl['peptide_outlier'].sum()), 0)
        self.assertEqual(int(pl['mainchain_outlier'].sum()), 0)

    def test_phi_psi_series(self):
        pps = self.ga.phi_psi_series()
        self.assertEqual(len(pps), 170)
        # N-terminal residue (PRO) has no phi and no previous residue
        self.assertIsNone(pps[0]['phi'])
        self.assertEqual(pps[0]['resname'], 'PRO')
        self.assertIsNone(pps[0]['prev_resname'])
        self.assertAlmostEqual(pps[1]['phi'], -70.19, places=2)

    def test_rama_molprobity(self):
        rm = self.ga.rama_molprobity()
        self.assertEqual(rm['analysed'], 168)   # 170 - 2 chain ends
        self.assertAlmostEqual(rm['favoured'], 0.917, places=3)
        self.assertAlmostEqual(rm['allowed'], 0.077, places=3)
        self.assertAlmostEqual(rm['outlier'], 0.006, places=3)


@ddt
class TestRamaGrid(TestCase):

    def setUp(self):
        # synthetic Top8000-style density: one dense bin, two sparse bins
        df = pd.DataFrame({
            'resname': ['ALA'] * 2500 + ['ALA', 'ALA'],
            'phi': [-62.0] * 2500 + [60.0, -135.0],
            'psi': [-45.0] * 2500 + [40.0, 140.0],
        })
        self.grid = RamaGrid.fit(df)

    @data(
        [None, None],
        [float('inf'), None],
        [-180.0, 0],
        [-0.1, 9],
        [179.9, 19],
    )
    @unpack
    def test_bin_index(self, deg, expect):
        self.assertEqual(RamaGrid._bin_index(deg), expect)

    @data(
        [['GLY'], 'GLY'],
        [['PRO'], 'PRO'],
        [['ALA'], 'general'],
        [['ALA', 'PRO'], 'pre-PRO'],
    )
    @unpack
    def test_rama_grid_key(self, args, expect):
        self.assertEqual(RamaGrid.rama_grid_key(*args), expect)

    def test_fit_all_grid_types_present(self):
        self.assertEqual(set(self.grid.grids.keys()),
                         set(RamaGrid.GRID_TYPES))
        for arr in self.grid.grids.values():
            self.assertEqual(arr.shape, (20, 20))

    def test_fit_cuts(self):
        # general: 2502 points -> favoured cut = the dense-bin count (2500),
        # allowed cut = the sparse-bin count (1); the empty bins (0) stay
        # below both cuts
        self.assertEqual(self.grid.cuts['general'], (2500, 1))
        # GLY grid empty -> no cuts
        self.assertEqual(self.grid.cuts['GLY'], (None, None))

    @data(
        [-62.0, -45.0, 'ALA', None, 'favoured'],    # dense bin
        [60.0, 40.0, 'ALA', None, 'allowed'],       # sparse bins
        [-135.0, 140.0, 'ALA', None, 'allowed'],
        [150.0, -150.0, 'ALA', None, 'outlier'],    # empty bin
        [None, None, 'ALA', None, 'NA'],
    )
    @unpack
    def test_classify(self, phi, psi, resname, prev, expect):
        self.assertEqual(self.grid.classify(phi, psi, resname, prev), expect)

    def test_classify_pre_pro(self):
        # ALA after PRO maps onto the (empty) pre-PRO grid
        self.assertEqual(self.grid.classify(-62.0, -45.0, 'ALA', 'PRO'),
                         'outlier')

    @skipUnless(os.path.exists(RAMA_NPZ), 'Rama8000 grid data not available')
    def test_load_npz(self):
        grid = RamaGrid.load_npz(RAMA_NPZ)
        self.assertIn('general', grid.grids)
        self.assertIsNotNone(grid.cuts['general'][0])


@ddt
class TestGeometry(TestCase):

    def setUp(self):
        self.geo = Geometry(PDB_6BSX)
        self.info = self.geo.run('A')

    def test_bond_length_angle_z(self):
        z = self.geo.bond_length_angle_z(self.geo.get_chain('A'))
        self.assertEqual(list(z.columns),
            ['res_pos', 'resname', 'kind', 'atoms', 'observed',
             'ideal', 'sd', 'z', 'outlier'])
        self.assertEqual(set(z['kind']), {'bond', 'angle'})
        self.assertFalse(z['outlier'].any())

    def test_chain_rmsz(self):
        rmsz = self.geo.chain_rmsz(self.geo.get_chain('A'))
        self.assertAlmostEqual(rmsz['bond_rmsz'], 0.6415, places=3)
        self.assertAlmostEqual(rmsz['angle_rmsz'], 0.7344, places=3)
        self.assertEqual(rmsz['n_bond_outliers'], 0)
        self.assertEqual(rmsz['n_angle_outliers'], 0)

    def test_run_summary(self):
        self.assertEqual(self.info['total_residues'], 170)
        self.assertEqual(self.info['num_missing'], 0)
        self.assertAlmostEqual(self.info['ca_distance']['mean'], 3.805,
                               places=3)
        self.assertAlmostEqual(self.info['peptide_distance']['max'], 1.361,
                               places=3)
        self.assertEqual(self.info['chain_continuity']['num_fragments'], 1)
        self.assertEqual(self.info['chain_continuity']['num_breaks'], 0)
        self.assertEqual(self.info['num_clashes'], 0)
        self.assertEqual(self.info['num_severe_clashes'], 0)
        self.assertAlmostEqual(self.info['radius_of_gyration'], 14.8,
                               places=1)
        self.assertAlmostEqual(self.info['omega']['cis'], 0.0)
        # raw per-metric data nested under 'details'
        self.assertIn('details', self.info)
        self.assertIn('missing_residues', self.info['details'])
        self.assertIn('chain_continuity', self.info['details'])

    def test_run_unknown_chain(self):
        self.assertEqual(self.geo.run('Z'), {})

    def test_check_pass(self):
        self.assertEqual(self.geo.check(self.info), (True, 'pass'))

    def test_check_strained_model_passes(self):
        # 2OX8: 2 clashes but none severe and a low clash score
        info = Geometry(PDB_2OX8).run('A')
        self.assertEqual(info['num_clashes'], 2)
        self.assertEqual(info['num_severe_clashes'], 0)
        self.assertEqual(Geometry(PDB_2OX8).check(info), (True, 'pass'))

    @data(
        # empty summary
        [{}, 'null_geometry'],
        # missing backbone atoms are not tolerated
        [{'num_missing': 2}, 'missing_residues'],
        # chain breaks
        [{'chain_continuity': {'num_breaks': 1, 'num_fragments': 2}},
         'chain_breaks'],
        # CA-CA stretched
        [{'ca_distance': {'max': 5.5}}, 'ca_distance'],
        # peptide bond broken / strained
        [{'peptide_distance': {'max': 5.5, 'mean': 1.4}}, 'peptide_distance'],
        [{'peptide_distance': {'max': 1.4, 'mean': 3.5}}, 'peptide_distance'],
        # too many unresolved torsions
        [{'torsion_angles': {'null': 0.2}}, 'torsion_angles'],
        # severe clash (< 2.0 A) vs plain clash score
        [{'num_severe_clashes': 1}, 'severe_clashes'],
        [{'clash_score': 31.0}, 'clashes'],
    )
    @unpack
    def test_check_gates(self, geo_info, expect_reason):
        passed, reason = self.geo.check(geo_info)
        self.assertFalse(passed)
        self.assertEqual(reason, expect_reason)

    def test_check_backbone_only_skips_torsions(self):
        backbone = Geometry(PDB_6BSX, is_backbone=True)
        passed, reason = backbone.check({'torsion_angles': {'null': 0.9}})
        self.assertTrue(passed)
        self.assertEqual(reason, 'pass')
