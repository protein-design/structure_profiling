'''
unit tests for chain-continuity analysis (qc/breaks.py, Breaks):
missing backbone atoms, chain breaks / fragments, strained peptides,
CA-CA / peptide C-N distances and the run() / check() summary + gates
- clean model: tests/data/PDB/6BSX_0_A_0.pdb
    (chain A: 170 residues; every peptide intact, one fragment -> pass)
- synthetic defects (all shipped PDB test data is clean, so they are
  built on the fly with the PDB-editing helpers in tests/helper.py):
    - broken:   C-terminal half (resseq >= 87) rigidly translated 20 A
                -> one break, two fragments ('chain_breaks' gate)
    - strained: peptide 86-87 stretched 0.5 A along C->N
                (1.5 A < d_C_N = 1.84 A <= 2.0 A, not broken)
    - missing:  O of residue 87 dropped -> num_missing = 1
    - gap:      residue 87 deleted -> numbering-gap break
                (resnum_gap = 2, 'chain_breaks' gate)
- pytest -s tests/qc/test_breaks.py
'''
import os
import tempfile

import numpy as np

from ..helper import *
from src.structure_profiling import Pdb, Breaks
from src.structure_profiling.geometry.geometry_distance import GeometryDistance

PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')

# residue where the synthetic defects are introduced
NUDGE_RES = 87
PAIR_COLUMNS = ['prev_res', 'prev_name', 'next_res', 'next_name',
                'resnum_gap', 'd_C_N', 'd_CA_CA']


class TestBreaks(TestCase):

    def setUp(self):
        self.chain = Pdb(PDB_6BSX).get_chain('A')
        self.breaks = Breaks(self.chain)
        self.lines = open(PDB_6BSX).readlines()
        # temp dir for the synthetic defect models
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def defect_pdb(self, name, lines):
        '''write edited PDB lines into the temp dir, return the path.'''
        return write_pdb_lines(lines, os.path.join(self.tmp.name, name))

    def broken_chain(self):
        '''C-terminal half rigidly translated 20 A: one break, two fragments.'''
        path = self.defect_pdb('broken.pdb',
            translate_residues(self.lines, NUDGE_RES, (20.0, 0.0, 0.0)))
        return Pdb(path).get_chain('A')

    def strained_chain(self):
        '''peptide 86-87 stretched into the strain window, not broken.'''
        path = self.defect_pdb('strained.pdb',
            stretch_peptide(self.lines, NUDGE_RES, extra=0.5))
        return Pdb(path).get_chain('A')

    def missing_chain(self):
        '''O of residue 87 dropped: one residue with missing backbone atoms.'''
        path = self.defect_pdb('missing.pdb',
            strip_atom(self.lines, NUDGE_RES, 'O'))
        return Pdb(path).get_chain('A')

    def gap_chain(self):
        '''residue 87 deleted: numbering-gap break (resnum_gap = 2).'''
        path = self.defect_pdb('gap.pdb', delete_residue(self.lines, NUDGE_RES))
        return Pdb(path).get_chain('A')

    # ------------------------------------------------------------------
    # selection + distances
    # ------------------------------------------------------------------

    def test_standard_residues(self):
        self.assertEqual(len(self.breaks.standard_residues()), 170)

    def test_missing_residues(self):
        # clean model: every residue has N, CA, C, O
        self.assertEqual(self.breaks.missing_residues(), [])
        missing = Breaks(self.missing_chain()).missing_residues()
        self.assertEqual(missing, [(87, {'O'})])

    def test_ca_and_peptide_distances(self):
        # ideal geometry of the clean model
        ca = self.breaks.ca_distances()
        pep = self.breaks.peptide_distances()
        self.assertEqual(len(ca), 169)
        self.assertEqual(len(pep), 169)
        self.assertAlmostEqual(float(np.mean(ca)), 3.805, places=3)

    # ------------------------------------------------------------------
    # chain continuity
    # ------------------------------------------------------------------

    def test_chain_continuity_clean(self):
        cont = self.breaks.chain_continuity()
        self.assertEqual(cont['num_residues'], 170)
        self.assertEqual(cont['num_fragments'], 1)
        self.assertEqual(cont['num_breaks'], 0)
        self.assertEqual(cont['num_strained'], 0)
        # column-complete, empty frames
        self.assertTrue(cont['breaks'].empty)
        self.assertEqual(list(cont['breaks'].columns),
                         PAIR_COLUMNS + ['fragment_after'])
        self.assertTrue(cont['strained'].empty)
        self.assertEqual(list(cont['strained'].columns), PAIR_COLUMNS)

    def test_chain_continuity_broken(self):
        cont = Breaks(self.broken_chain()).chain_continuity()
        self.assertEqual(cont['num_breaks'], 1)
        self.assertEqual(cont['num_fragments'], 2)
        row = cont['breaks'].iloc[0]
        self.assertEqual(int(row['prev_res']), 86)
        self.assertEqual(row['prev_name'], 'SER')
        self.assertEqual(int(row['next_res']), 87)
        self.assertEqual(row['next_name'], 'VAL')
        self.assertEqual(int(row['resnum_gap']), 1)
        self.assertAlmostEqual(float(row['d_C_N']), 18.93, places=2)
        self.assertAlmostEqual(float(row['d_CA_CA']), 16.22, places=2)
        self.assertEqual(int(row['fragment_after']), 1)

    def test_chain_continuity_numbering_gap(self):
        cont = Breaks(self.gap_chain()).chain_continuity()
        self.assertEqual(cont['num_residues'], 169)
        self.assertEqual(cont['num_breaks'], 1)
        self.assertEqual(cont['num_fragments'], 2)
        row = cont['breaks'].iloc[0]
        # residue 87 deleted: 86 -> 88 skips one number and the peptide
        # is unmeasurably far (d_C_N ~ 2.99 A)
        self.assertEqual((int(row['prev_res']), int(row['next_res'])), (86, 88))
        self.assertEqual(int(row['resnum_gap']), 2)
        self.assertAlmostEqual(float(row['d_C_N']), 2.99, places=2)

    def test_chain_continuity_strained(self):
        cont = Breaks(self.strained_chain()).chain_continuity()
        # stretched but intact: not a break, no new fragment
        self.assertEqual(cont['num_breaks'], 0)
        self.assertEqual(cont['num_fragments'], 1)
        self.assertEqual(cont['num_strained'], 1)
        row = cont['strained'].iloc[0]
        self.assertEqual((int(row['prev_res']), int(row['next_res'])), (86, 87))
        # Breaks.cn_strain < d_C_N <= Breaks.cn_break
        self.assertGreater(float(row['d_C_N']), Breaks.cn_strain)
        self.assertLessEqual(float(row['d_C_N']), Breaks.cn_break)
        self.assertAlmostEqual(float(row['d_C_N']), 1.84, places=2)


    # ------------------------------------------------------------------
    # run() summary + check() gates
    # ------------------------------------------------------------------

    def test_run_summary(self):
        info = self.breaks.run(keep_details=False)
        self.assertEqual(set(info), {'n_residues', 'num_missing',
            'num_breaks', 'num_fragments', 'num_strained', 'ca_max_dist',
            'pep_max_dist'})
        self.assertEqual(info['n_residues'], 170)
        self.assertEqual(info['num_missing'], 0)
        self.assertEqual(info['num_breaks'], 0)
        self.assertEqual(info['num_fragments'], 1)
        self.assertEqual(info['num_strained'], 0)
        self.assertEqual(info['ca_max_dist'], 3.88)
        self.assertEqual(info['pep_max_dist'], 1.36)
        # details hold the raw tables
        detailed = self.breaks.run()
        self.assertIn('missing_residues', detailed['details'])
        self.assertIn('breaks', detailed['details'])
        self.assertIn('strained', detailed['details'])

    def test_run_broken(self):
        info = Breaks(self.broken_chain()).run(keep_details=False)
        self.assertEqual(info['num_breaks'], 1)
        self.assertEqual(info['num_fragments'], 2)
        self.assertEqual(info['ca_max_dist'], 16.22)
        self.assertEqual(info['pep_max_dist'], 18.93)

    def test_check(self):
        # clean model passes all gates
        self.assertEqual(self.breaks.check(), (True, 'pass'))
        # real break and numbering gap: 'chain_breaks'
        self.assertEqual(Breaks(self.broken_chain()).check(),
                         (False, 'chain_breaks'))
        self.assertEqual(Breaks(self.gap_chain()).check(),
                         (False, 'chain_breaks'))
        # missing backbone atoms are not tolerated
        self.assertEqual(Breaks(self.missing_chain()).check(),
                         (False, 'missing_residues'))
        # one strained peptide is within tolerance (gate: num_strained > 3)
        self.assertEqual(Breaks(self.strained_chain()).check(),
                         (True, 'pass'))
        # synthetic gates: fragments, strain, distances, null summary
        self.assertEqual(self.breaks.check({'num_fragments': 2}),
                         (False, 'chain_fragments'))
        self.assertEqual(self.breaks.check({'num_strained': 4}),
                         (False, 'chain_strain'))
        self.assertEqual(self.breaks.check({'ca_max_dist': 6.0}),
                         (False, 'ca_distance'))
        self.assertEqual(self.breaks.check({'pep_max_dist': 6.0}),
                         (False, 'peptide_distance'))
        self.assertEqual(self.breaks.check({}), (False, 'null_breaks'))

    # ------------------------------------------------------------------
    # backward-compatible delegates on GeometryDistance
    # ------------------------------------------------------------------

    def test_geometry_delegates_agree(self):
        # GeometryDistance continuity APIs must stay consistent with Breaks
        gd = GeometryDistance(self.chain)
        self.assertEqual(gd.missing_residues(),
                         self.breaks.missing_residues())
        np.testing.assert_allclose(gd.ca_distances(),
                                   self.breaks.ca_distances())
        np.testing.assert_allclose(gd.peptide_distances(),
                                   self.breaks.peptide_distances())
        self.assertTrue(gd.chain_continuity()['breaks'].equals(
            self.breaks.chain_continuity()['breaks']))
        # same agreement on the broken model
        gd_b = GeometryDistance(self.broken_chain())
        brk_b = Breaks(self.broken_chain())
        self.assertTrue(gd_b.chain_continuity()['breaks'].equals(
            brk_b.chain_continuity()['breaks']))
        self.assertEqual(gd_b.missing_residues(), brk_b.missing_residues())

