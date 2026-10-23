'''
unit tests for chain-level torsion / angle QC (qc/torsions.py, Torsions):
bond-angle Z-scores vs Engh & Huber ideals, phi/psi/omega torsions with
cis/trans conformation, wwPDB peptide / main-chain planarity, Ramachandran
(single-box + MolProbity style + RamaGrid) and the run() / check()
summary + gates
- clean model: tests/data/PDB/6BSX_0_A_0.pdb
    (chain A: 170 residues, numbering 2..171; every backbone angle within
    |Z| <= 5, planarity clean, one MolProbity Rama outlier (ASP 65)
    -> pass)
- synthetic defects (all shipped PDB test data is clean, so they are
  built on the fly with the PDB-editing helpers in tests/helper.py):
    - rotated:  residues with resseq >= 10 rigidly rotated 180 deg about
                their own CA -> intra-residue geometry unchanged, but the
                160 peptide bridges downstream scramble: CA-C-N / C-N-CA
                angles leave the Engh & Huber window (n_z_outliers = 320
                = 2 per peptide), omega and phi / psi at the junctions
                become garbage ('angle_rmsz' gate)
    - missing:  O of residue 10 dropped -> peptide 10-11 is not
                measurable (n_peptides = 168) but nothing else changes
                (no gate; phi / psi only need N, CA, C)
    - angles:   N of residue 8 pulled 4 A from its CA -> exactly the two
                bridge bond angles at that N collapse (28.3 / 31.6 deg,
                |Z| = 44 / 35) -> 'angle_rmsz' gate
- the GeometryAngle delegates on the moved methods must agree with
  Torsions (backward compatibility)
- pytest -s tests/qc/test_torsions.py
'''
import math
import os
import tempfile

import numpy as np
import pandas as pd

from ..helper import *
from src.structure_profiling import Pdb, Torsions
from src.structure_profiling.qc.geometry_angle import GeometryAngle
from src.structure_profiling.qc.rama_grid import RamaGrid

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')

# residues where the synthetic defects are introduced
NUDGE_RES = 10     # rotated / missing-O defects
ANGLE_RES = 8      # pulled-N bond-angle defect


class TestTorsions(TestCase):

    def setUp(self):
        self.chain = Pdb(PDB_6BSX).get_chain('A')
        self.torsions = Torsions(self.chain)
        self.lines = open(PDB_6BSX).readlines()
        # temp dir for the synthetic defect models
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def defect_pdb(self, name, lines):
        '''write edited PDB lines into the temp dir, return the path.'''
        return write_pdb_lines(lines, os.path.join(self.tmp.name, name))

    def rotated_chain(self):
        '''residues with resseq >= 10 rigidly rotated 180 deg about
        their own CA atoms.'''
        path = self.defect_pdb('rotated.pdb',
            rotate_residues(self.lines, NUDGE_RES, 180.0))
        return Pdb(path).get_chain('A')

    def missing_chain(self):
        '''O of residue 10 dropped: peptide 10-11 not measurable.'''
        path = self.defect_pdb('missing.pdb',
            strip_atom(self.lines, NUDGE_RES, 'O'))
        return Pdb(path).get_chain('A')

    def angle_chain(self):
        '''N of residue 8 pulled radially to 4 A from its CA: the two
        bridge bond angles at that N collapse, the rest of the chain is
        untouched.'''
        res = {r.id[1]: r for r in self.torsions.standard_residues()}
        r = res[ANGLE_RES]
        ca, n = r['CA'].coord, r['N'].coord
        target = ca + 4.0 * (n - ca) / np.linalg.norm(n - ca)
        path = self.defect_pdb('angles.pdb',
            move_atom(self.lines, ANGLE_RES, 'N', target))
        return Pdb(path).get_chain('A')

    # ------------------------------------------------------------------
    # point-based angle utilities
    # ------------------------------------------------------------------

    def test_dihedral_and_bond_angle(self):
        p = lambda *c: np.array(c, dtype=float)
        self.assertAlmostEqual(Torsions.dihedral_angle(
            p(0, 1, 1), p(0, 0, 0), p(1, 0, 0), p(1, 1, 1)), 0.0)
        self.assertAlmostEqual(Torsions.dihedral_angle(
            p(0, 1, 1), p(0, 0, 0), p(1, 0, 0), p(1, -1, -1)), 180.0)
        self.assertAlmostEqual(Torsions.bond_angle(
            p(1, 0, 0), p(0, 0, 0), p(0, 1, 0)), 90.0)

    def test_evaluate_ramachandran(self):
        self.assertEqual(Torsions.evaluate_ramachandran(-60, -45),
                         'right-handed-helix')
        self.assertEqual(Torsions.evaluate_ramachandran(-135, 140),
                         'beta-sheet')
        self.assertEqual(Torsions.evaluate_ramachandran(None, None), 'null')

    def test_rama_class(self):
        self.assertEqual(Torsions.rama_class(-62, -41, 'ALA'), 'favoured')
        self.assertEqual(Torsions.rama_class(-100, 0, 'PRO'), 'allowed')
        self.assertEqual(Torsions.rama_class(0, 0, 'PRO'), 'outlier')
        self.assertEqual(Torsions.rama_class(None, None, 'ALA'), 'NA')

    def test_ideal_angles_shared(self):
        # canonical values live in GeometryAngle, Torsions re-exports them
        self.assertIs(Torsions.IDEAL_ANGLES, GeometryAngle.IDEAL_ANGLES)
        self.assertIs(Torsions.Z_OUTLIER_CUT, GeometryAngle.Z_OUTLIER_CUT)

    # ------------------------------------------------------------------
    # backbone bond angles
    # ------------------------------------------------------------------

    def test_bond_angles(self):
        ba = self.torsions.bond_angles()
        self.assertEqual(len(ba), 170)
        self.assertAlmostEqual(ba['angle_N_CA_C'].mean(), 111.06, places=2)
        self.assertEqual(int(ba['angle_CA_C_N'].notna().sum()), 169)

    def test_bond_angle_z_clean(self):
        az = self.torsions.bond_angle_z()
        # 170 N-CA-C + 2 * 169 bridge angles
        self.assertEqual(len(az), 508)
        self.assertEqual(az['atoms'].value_counts().to_dict(),
                         {'N-CA-C': 170, 'CA-C-N': 169, 'C-N-CA': 169})
        self.assertEqual(int(az['outlier'].sum()), 0)

    def test_bond_angle_z_pulled_n(self):
        # the moved N collapses exactly its two bridge angles; the
        # intra-residue N-CA-C angle is untouched (radial move)
        az = Torsions(self.angle_chain()).bond_angle_z()
        bad = az.loc[az['outlier']]
        self.assertEqual(len(bad), 2)
        self.assertEqual(bad['res_pos'].tolist(), [7, 7])
        self.assertEqual(bad['atoms'].tolist(), ['CA-C-N', 'C-N-CA'])
        self.assertAlmostEqual(bad['observed'].iloc[0], 28.31, places=2)
        self.assertAlmostEqual(bad['z'].iloc[0], -43.94, places=2)
        self.assertAlmostEqual(bad['observed'].iloc[1], 31.57, places=2)
        self.assertAlmostEqual(bad['z'].iloc[1], -34.74, places=2)

    def test_bond_angle_z_rotated_preserves_intra(self):
        # rigid per-residue rotation keeps every N-CA-C angle in range
        # while the 160 rotated-side bridges (CA-C-N / C-N-CA) blow up
        az = Torsions(self.rotated_chain()).bond_angle_z()
        bad = az.loc[az['outlier']]
        self.assertEqual(len(bad), 320)
        self.assertEqual((bad['atoms'] == 'N-CA-C').sum(), 0)

    # ------------------------------------------------------------------
    # backbone torsions (phi / psi / omega)
    # ------------------------------------------------------------------

    def test_torsion_angles(self):
        ta = self.torsions.torsion_angles()
        self.assertEqual(len(ta), 170)
        # numbering 2..171: N-terminal phi and C-terminal psi / omega
        # are not measurable
        self.assertEqual(ta['res_pos'].iloc[0], 2)
        self.assertEqual(ta['residue'].iloc[0], 'PRO')
        self.assertTrue(math.isnan(ta['phi'].iloc[0]))
        self.assertAlmostEqual(ta['omega'].iloc[0], -178.99, places=2)
        self.assertAlmostEqual(ta['phi'].iloc[1], -70.19, places=2)
        self.assertTrue(math.isnan(ta['psi'].iloc[-1]))
        self.assertTrue(math.isnan(ta['omega'].iloc[-1]))
        self.assertEqual(ta['conformation'].value_counts().to_dict(),
                         {'trans': 169, 'null': 1})
        # single-box Ramachandran labels (evaluate_ramachandran)
        self.assertEqual(ta['ramachandran'].value_counts().to_dict(),
                         {'right-handed-helix': 80, 'beta-sheet': 68,
                          'allowed': 6, 'bad': 9, 'left-handed-helix': 5,
                          'null': 2})

    def test_torsion_angles_rotated(self):
        ta = Torsions(self.rotated_chain()).torsion_angles()
        # peptides up to 8-9 keep their omega, the first rotated-side
        # peptide (9-10) leaves the trans window
        self.assertAlmostEqual(ta['omega'].iloc[6], -179.91, places=2)
        self.assertAlmostEqual(ta['omega'].iloc[7], 135.44, places=2)
        # rigid rotation keeps N / CA / C in place: no null phi / psi
        self.assertEqual(ta['conformation'].value_counts().to_dict(),
                         {'trans': 169, 'null': 1})

    # ------------------------------------------------------------------
    # wwPDB peptide / main-chain planarity
    # ------------------------------------------------------------------

    def test_peptide_planarity(self):
        pl = self.torsions.peptide_planarity()
        self.assertEqual(len(pl), 169)
        self.assertEqual(int(pl['peptide_outlier'].sum()), 0)
        self.assertEqual(int(pl['mainchain_outlier'].sum()), 0)
        self.assertAlmostEqual(pl['omega_dev'].iloc[0], 1.01, places=2)

    def test_peptide_planarity_missing_o(self):
        # peptide 10-11 (CA,C,O of 10 + N,CA of 11) is skipped without
        # O of 10; the remaining peptides stay measurable and clean
        pl = Torsions(self.missing_chain()).peptide_planarity()
        self.assertEqual(len(pl), 168)
        self.assertNotIn(NUDGE_RES, pl['res_pos'].values)
        self.assertEqual(int(pl['peptide_outlier'].sum()), 0)
        self.assertEqual(int(pl['mainchain_outlier'].sum()), 0)

    # ------------------------------------------------------------------
    # Ramachandran: single-box, MolProbity style, RamaGrid
    # ------------------------------------------------------------------

    def test_phi_psi_series(self):
        pps = self.torsions.phi_psi_series()
        self.assertEqual(len(pps), 170)
        # N-terminal residue (PRO 2) has no phi and no previous residue
        self.assertIsNone(pps[0]['phi'])
        self.assertEqual(pps[0]['resname'], 'PRO')
        self.assertIsNone(pps[0]['prev_resname'])
        self.assertAlmostEqual(pps[1]['phi'], -70.19, places=2)

    def test_rama_molprobity(self):
        rm = self.torsions.rama_molprobity()
        self.assertEqual(rm['analysed'], 168)   # 170 - phi/psi chain ends
        self.assertAlmostEqual(rm['favoured'], 0.917, places=3)
        self.assertAlmostEqual(rm['allowed'], 0.077, places=3)
        self.assertAlmostEqual(rm['outlier'], 0.006, places=3)
        # the single outlier: ASP 65 stretched beta (phi -143, psi -167)
        self.assertEqual([(o[0], o[1]) for o in rm['outliers']],
                         [(65, 'ASP')])
        self.assertAlmostEqual(rm['outliers'][0][2], -143.38, places=2)
        self.assertAlmostEqual(rm['outliers'][0][3], -166.88, places=2)

    def test_rama_molprobity_rotated(self):
        # half-rotated chain: outliers explode, favoured collapses
        rm = Torsions(self.rotated_chain()).rama_molprobity()
        self.assertAlmostEqual(rm['outlier'], 0.476, places=3)
        self.assertAlmostEqual(rm['favoured'], 0.131, places=3)
        self.assertAlmostEqual(rm['allowed'], 0.393, places=3)

    def test_rama_grid_summary(self):
        # synthetic Top8000-style density: one dense alpha-helix bin
        df = pd.DataFrame({
            'resname': ['ALA'] * 2500,
            'phi': [-62.0] * 2500,
            'psi': [-45.0] * 2500,
        })
        gs = self.torsions.rama_grid_summary(RamaGrid.fit(df))
        self.assertEqual(set(gs.keys()),
                         {'analysed', 'favoured', 'allowed', 'outlier',
                          'outliers'})
        # the two chain ends (phi / psi None -> NA) are not analysed
        self.assertEqual(gs['analysed'], 168)
        # one dense bin only: residues in it are favoured, every other
        # measured residue falls on an empty bin -> outlier, nothing
        # lands in between (allowed = 0)
        self.assertEqual(gs['allowed'], 0.0)
        self.assertGreater(gs['favoured'], 0.0)
        self.assertAlmostEqual(gs['favoured'] + gs['outlier'], 1.0,
                               places=2)
        # the outliers carry residue context; chain ends are NA
        self.assertEqual(len(gs['outliers']),
                         round(gs['outlier'] * 168))
        self.assertNotIn(2, [o[0] for o in gs['outliers']])
        self.assertNotIn(171, [o[0] for o in gs['outliers']])

    # ------------------------------------------------------------------
    # summary + gates (chain_torsions table)
    # ------------------------------------------------------------------

    def test_run_clean(self):
        info = self.torsions.run(keep_details=False)
        self.assertEqual(set(info.keys()), {   # the chain_torsions columns
            'n_residues', 'n_angles', 'angle_rmsz', 'n_z_outliers',
            'z_outlier_frac', 'n_torsions', 'torsion_null_frac',
            'frac_trans', 'frac_cis', 'n_cis_nonpro', 'rama_analysed',
            'rama_favoured_frac', 'rama_allowed_frac', 'rama_outlier_frac',
            'n_peptides', 'peptide_outlier_frac', 'mainchain_outlier_frac'})
        self.assertEqual(info['n_residues'], 170)
        self.assertEqual(info['n_angles'], 508)
        self.assertAlmostEqual(info['angle_rmsz'], 0.734, places=3)
        self.assertEqual(info['n_z_outliers'], 0)
        self.assertEqual(info['z_outlier_frac'], 0.0)
        self.assertEqual(info['n_torsions'], 170)
        self.assertAlmostEqual(info['torsion_null_frac'], 0.012, places=3)
        self.assertEqual(info['frac_trans'], 1.0)
        self.assertEqual(info['frac_cis'], 0.0)
        self.assertEqual(info['n_cis_nonpro'], 0)
        self.assertEqual(info['rama_analysed'], 168)
        self.assertAlmostEqual(info['rama_favoured_frac'], 0.917, places=3)
        self.assertAlmostEqual(info['rama_allowed_frac'], 0.077, places=3)
        self.assertAlmostEqual(info['rama_outlier_frac'], 0.006, places=3)
        self.assertEqual(info['n_peptides'], 169)
        self.assertEqual(info['peptide_outlier_frac'], 0.0)
        self.assertEqual(info['mainchain_outlier_frac'], 0.0)
        # details hold the raw tables
        detailed = self.torsions.run()
        self.assertEqual(set(detailed['details'].keys()),
                         {'bond_angle_z', 'torsion_angles',
                          'peptide_planarity', 'rama_outliers'})

    def test_run_missing_o(self):
        # a missing O only shrinks the peptide table, every other
        # statistic is unchanged
        info = Torsions(self.missing_chain()).run(keep_details=False)
        self.assertEqual(info['n_peptides'], 168)
        self.assertEqual(info['peptide_outlier_frac'], 0.0)
        self.assertEqual(info['n_torsions'], 170)
        self.assertAlmostEqual(info['angle_rmsz'], 0.734, places=3)

    def test_run_rotated(self):
        info = Torsions(self.rotated_chain()).run(keep_details=False)
        self.assertAlmostEqual(info['angle_rmsz'], 29.588, places=3)
        self.assertEqual(info['n_z_outliers'], 320)
        self.assertAlmostEqual(info['z_outlier_frac'], 0.63, places=3)
        self.assertAlmostEqual(info['peptide_outlier_frac'], 0.065, places=3)
        self.assertAlmostEqual(info['mainchain_outlier_frac'], 0.562,
                               places=3)
        self.assertAlmostEqual(info['rama_outlier_frac'], 0.476, places=3)
        # torsions themselves stay measurable: rigid rotation keeps
        # N / CA / C of every residue in place
        self.assertEqual(info['n_torsions'], 170)
        self.assertAlmostEqual(info['torsion_null_frac'], 0.012, places=3)

    def test_run_bad_angle(self):
        info = Torsions(self.angle_chain()).run(keep_details=False)
        self.assertAlmostEqual(info['angle_rmsz'], 2.591, places=3)
        self.assertEqual(info['n_z_outliers'], 2)
        self.assertAlmostEqual(info['z_outlier_frac'], 0.004, places=3)

    def test_check(self):
        # clean model passes all gates
        self.assertEqual(self.torsions.check(), (True, 'pass'))
        # 320 collapsed bridge angles and the pulled-N pair: angle rmsz
        self.assertEqual(Torsions(self.rotated_chain()).check(),
                         (False, 'angle_rmsz'))
        self.assertEqual(Torsions(self.angle_chain()).check(),
                         (False, 'angle_rmsz'))
        # a missing O only shrinks the peptide table (within tolerance)
        self.assertEqual(Torsions(self.missing_chain()).check(),
                         (True, 'pass'))
        # synthetic summary overrides for every gate
        self.assertEqual(self.torsions.check({'angle_rmsz': 2.0}),
                         (False, 'angle_rmsz'))
        self.assertEqual(self.torsions.check({'z_outlier_frac': 0.02}),
                         (False, 'z_outliers'))
        self.assertEqual(self.torsions.check({'torsion_null_frac': 0.1}),
                         (False, 'torsion_null'))
        self.assertEqual(self.torsions.check({'peptide_outlier_frac': 0.1}),
                         (False, 'peptide_planarity'))
        self.assertEqual(self.torsions.check(
            {'mainchain_outlier_frac': 0.1}), (False, 'mainchain_planarity'))
        self.assertEqual(self.torsions.check({'rama_outlier_frac': 0.1}),
                         (False, 'rama_outliers'))
        # null summary
        self.assertEqual(self.torsions.check({}), (False, 'null_torsions'))

    def test_gates_and_cutoffs(self):
        self.assertEqual(Torsions.GATES, {
            'angle_rmsz': 1.5,
            'z_outlier_frac': 0.01,
            'torsion_null_frac': 0.05,
            'peptide_outlier_frac': 0.05,
            'mainchain_outlier_frac': 0.05,
            'rama_outlier_frac': 0.05})
        self.assertEqual(Torsions.Z_OUTLIER_CUT, 5.0)
        self.assertEqual(Torsions.CIS_CUT, 90.0)
        self.assertEqual(Torsions.OMEGA_DEV_CUT, 30.0)
        self.assertEqual(Torsions.MAINCHAIN_OOP_CUT, 10.0)

    # ------------------------------------------------------------------
    # backward-compatible delegates on GeometryAngle
    # ------------------------------------------------------------------

    def test_geometry_delegates_agree(self):
        # GeometryAngle APIs moved to Torsions must stay consistent
        ga = GeometryAngle(self.chain)
        for meth in ('bond_angles', 'torsion_angles', 'peptide_planarity',
                     'bond_angle_z'):
            self.assertTrue(getattr(ga, meth)().equals(
                getattr(self.torsions, meth)()), meth)
        self.assertEqual(ga.phi_psi_series(), self.torsions.phi_psi_series())
        self.assertEqual(ga.rama_molprobity(),
                         self.torsions.rama_molprobity())
        self.assertEqual(ga.evaluate_ramachandran(-60, -45),
                         Torsions.evaluate_ramachandran(-60, -45))
        self.assertEqual(ga.rama_class(-62, -41, 'ALA'),
                         Torsions.rama_class(-62, -41, 'ALA'))
        # same agreement on the rotated model
        ga_r = GeometryAngle(self.rotated_chain())
        to_r = Torsions(self.rotated_chain())
        self.assertTrue(ga_r.bond_angle_z().equals(to_r.bond_angle_z()))
        self.assertTrue(ga_r.torsion_angles().equals(
            to_r.torsion_angles()))
