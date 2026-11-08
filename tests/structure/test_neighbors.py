'''
test the Neighbors ligand / binding-site neighborhood search on a small
synthetic chain (exact coordinates, see tests.helper.synthesize_site_pdb)
plus the shipped 2OX8 chain PDB (no HETATM -> empty result).

at the default 5.0 A cutoff the synthetic ligand LIG 300 touches exactly
residues ALA 2 and ALA 102 (min CA-C1 4.5 A): neighbors
close in space but scattered over the sequence (n_fragments = 2).
- pytest -s tests/structure/test_neighbors.py
'''
import math
import os
import tempfile
from pathlib import Path
from unittest import TestCase, skipUnless

from src.structure_profiling import Pdb, Neighbors

from ..helper import DATA_DIR, synthesize_site_pdb

PDB_2OX8 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')

# exact expected contact distances (A) of the synthetic fixture
MIN_DIST_2 = math.sqrt(16.25)      # CB(ALA 2/102) - C1 ligand
MIN_DIST_SITE = math.sqrt(3.4)     # C(ALA 1/3) - N(ALA 2) active site


def neighbors_structure():
    '''load the synthetic site PDB (written to a temp dir).'''
    tmp = tempfile.TemporaryDirectory()
    pdb_file = synthesize_site_pdb(Path(tmp.name) / 'site.pdb')
    p = Pdb(pdb_file, None, False)
    tmp.cleanup()   # structure is fully parsed into memory
    return p.structure


def nb_fragments(df):
    '''Neighbors.fragments on a deduplicated copy of df.'''
    uniq = df.drop_duplicates(['chain_no', 'resseq'])
    return Neighbors(None).fragments(uniq)


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestNeighborsLigandSites(TestCase):
    '''site resolution: HETATM ligands, waters, explicit residues'''

    def test_ligand_sites_default(self):
        nb = Neighbors(neighbors_structure())
        sites = nb.ligand_sites()
        self.assertEqual([label for label, _ in sites], ['A/LIG300'])
        self.assertEqual([a.name for a in sites[0][1]], ['C1', 'N1'])

    def test_ligand_sites_include_waters(self):
        nb = Neighbors(neighbors_structure())
        sites = nb.ligand_sites(include_waters=True)
        self.assertEqual([label for label, _ in sites],
                         ['A/LIG300', 'A/HOH501'])

    def test_site_residues_activity_site(self):
        nb = Neighbors(neighbors_structure())
        sites = nb.site_residues([('A', 2)])
        self.assertEqual([label for label, _ in sites], ['A/ALA2'])
        self.assertEqual(len(sites[0][1]), 5)   # N, CA, C, O, CB

    def test_site_residues_missing_is_skipped(self):
        nb = Neighbors(neighbors_structure())
        self.assertEqual(nb.site_residues([('A', 999)]), [])

    def test_no_ligand_sites_without_hetatm(self):
        p = Pdb(PDB_2OX8, None, False)
        nb = Neighbors(p.structure)
        self.assertEqual(nb.ligand_sites(), [])
        self.assertEqual(nb.ligand_sites(include_waters=True), [])

    def test_chain_no_restricts_ligand_sites(self):
        # LIG300 / HOH501 sit in chain A: another chain has no sites
        nb = Neighbors(neighbors_structure(), chain_no='B')
        self.assertEqual(nb.ligand_sites(), [])
        self.assertEqual(nb.ligand_sites(include_waters=True), [])


@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestNeighborsSearch(TestCase):
    '''per-residue contact table of one site'''

    def setUp(self):
        self.nb = Neighbors(neighbors_structure())
        self.site = self.nb.ligand_sites()[0]

    def test_columns_and_order(self):
        df = self.nb.search(self.site[1])
        self.assertEqual(list(df.columns),
                         ['chain_no', 'resseq', 'icode', 'resname',
                          'min_dist', 'n_contacts', 'atom',
                          'site_resname', 'site_resseq', 'site_atom'])
        self.assertEqual(df['resseq'].tolist(), [2, 3, 102, 103])
        # the site residue (name + position) on every row
        self.assertEqual(df['site_resname'].unique().tolist(), ['LIG'])
        self.assertEqual(df['site_resseq'].unique().tolist(), [300])

    def test_exact_distances_and_contacts(self):
        df = self.nb.search(self.site[1])
        row = df.iloc[0]
        self.assertEqual(row['chain_no'], 'A')
        self.assertEqual(row['resname'], 'ALA')
        self.assertEqual(row['min_dist'], round(MIN_DIST_2, 3))
        self.assertEqual(row['atom'], 'CB')
        self.assertEqual(row['site_resname'], 'LIG')
        self.assertEqual(row['site_resseq'], 300)
        self.assertEqual(row['site_atom'], 'C1')
        # N, CA, C, CB within 5.0 A (O is at 5.669 A)
        self.assertEqual(row['n_contacts'], 4)

    def test_cutoff_widens_the_pocket(self):
        df = self.nb.search(self.site[1], cutoff=6.0)
        self.assertEqual(df['resseq'].tolist(), [1, 2, 3, 101, 102, 103])
        self.assertEqual(nb_fragments(df), 2)

    def test_chain_no_restricts_candidates(self):
        nb = Neighbors(self.nb.structure, chain_no='B')
        self.assertTrue(nb.search(self.site[1]).empty)

    def test_search_without_site_or_protein(self):
        self.assertTrue(self.nb.search([]).empty)
        empty = Neighbors(Pdb(PDB_2OX8, None, False).structure,
                          chain_no='Z')
        self.assertTrue(empty.search(self.site[1]).empty)

@skipUnless(os.path.isfile(PDB_2OX8), 'requires tests/data/PDB/2OX8_A.pdb')
class TestNeighborsRun(TestCase):
    '''summary statistics exported to chain_ligand'''

    def setUp(self):
        self.nb = Neighbors(neighbors_structure())

    def test_run_summary(self):
        info = self.nb.run(keep_details=False)
        self.assertEqual(info['n_sites'], 1)
        self.assertEqual(info['ligand_names'], 'LIG')
        self.assertEqual(info['n_neighbor_residues'], 4)
        self.assertEqual(info['n_fragments'], 2)
        self.assertAlmostEqual(info['min_dist'], MIN_DIST_2, places=3)
        self.assertAlmostEqual(info['mean_dist'], 4.45, places=2)
        self.assertEqual(info['min_dist_residues'],
                         'A/ALA2:4.031,A/ALA102:4.031')
        self.assertEqual(info['neighbor_residues'],
                         'A/ALA2:4.031,A/ALA3:4.867,'
                         'A/ALA102:4.031,A/ALA103:4.867')
        self.assertNotIn('neighbors', info)

    def test_run_details(self):
        info = self.nb.run()
        self.assertEqual(len(info['neighbors']), 4)
        self.assertIn('site', info['neighbors'].columns)

    def test_run_include_waters_adds_site(self):
        info = self.nb.run(include_waters=True, keep_details=False)
        self.assertEqual(info['n_sites'], 2)
        self.assertEqual(info['ligand_names'], 'LIG,HOH')
        # the water (x = 7.6) touches residues 3 and 103 instead
        self.assertEqual(info['n_neighbor_residues'], 4)

    def test_run_activity_site(self):
        info = self.nb.run(site_residues=[('A', 2)], keep_details=False)
        self.assertEqual(info['n_sites'], 1)
        self.assertEqual(info['ligand_names'], 'ALA')
        self.assertEqual(info['n_neighbor_residues'], 2)
        self.assertEqual(info['n_fragments'], 2)    # residues 1 and 3
        self.assertAlmostEqual(info['min_dist'], round(MIN_DIST_SITE, 3),
                               places=3)
        self.assertEqual(info['min_dist_residues'],
                         'A/ALA1:1.844,A/ALA3:1.844')

    def test_run_without_ligands(self):
        p = Pdb(PDB_2OX8, None, False)
        info = Neighbors(p.structure).run(keep_details=False)
        self.assertEqual(info['n_sites'], 0)
        self.assertEqual(info['n_neighbor_residues'], 0)
        self.assertEqual(info['n_fragments'], 0)
        self.assertIsNone(info['min_dist'])
        self.assertEqual(info['min_dist_residues'], '')
        self.assertEqual(info['neighbor_residues'], '')

    def test_run_without_structure(self):
        info = Neighbors(None).run(keep_details=False)
        self.assertEqual(info['n_sites'], 0)
        self.assertEqual(info['n_neighbor_residues'], 0)
