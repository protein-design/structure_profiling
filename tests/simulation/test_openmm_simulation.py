'''
integration tests for OpenmmSimulation on a small fixed AlphaFold model
- relies on OpenmmPdbfixer output (fixed.pdb with hydrogens)
- pytest -s tests/simulation/test_openmm_simulation.py
'''
import gzip
import os
import shutil
import tempfile
from pathlib import Path

import openmm.unit as openmm_unit

from ..helper import *
from src.structure_profiling.structure.openmm_pdbfixer import OpenmmPdbfixer
from src.structure_profiling.structure.openmm_simulation import OpenmmSimulation

PDB_LYSOZYME = os.path.join(
    DATA_DIR, 'lysozyme', 'AF-P61626-F1-model_v4.pdb')


@ddt
class TestOpenmmSimulation(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir = Path(tempfile.mkdtemp(prefix='bioomics_sim_'))
        # a fixed structure with hydrogens to simulate
        fixer = OpenmmPdbfixer(Path(PDB_LYSOZYME), {'outdir': cls.outdir})
        info = fixer.fast_fixer()
        assert info['status'] == 'succeed', info
        cls.fixed_pdb = str(cls.outdir / 'fixed.pdb')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir, ignore_errors=True)

    def test_init(self):
        obj = OpenmmSimulation(
            PDB_LYSOZYME, {'outdir': self.outdir / 'sim'})
        self.assertIsNone(obj.pdb)
        self.assertTrue(Path(obj.outdir).exists())

    def test_load_pdb(self):
        obj = OpenmmSimulation(
            PDB_LYSOZYME, {'outdir': self.outdir / 'sim'})
        obj.load_pdb()
        self.assertGreater(obj.pdb.topology.getNumAtoms(), 0)
        # the AF model has no CRYST1 record: a default box (>= 2x the
        # 1 nm nonbonded cutoff) must have been added
        self.assertEqual(obj.info['periodic_box'], 'default')
        box = obj.pdb.topology.getPeriodicBoxVectors()
        trace = (box[0][0] + box[1][1] + box[2][2]).value_in_unit(
            openmm_unit.nanometer)
        self.assertGreaterEqual(trace, 6.0)

    def test_load_pdb_gzipped(self):
        # regression: app.PDBFile() has no 'pdbfile' keyword and the .gz
        # check used Path.stem, which strips '.gz'
        gz_file = self.outdir / 'structure.pdb.gz'
        with open(PDB_LYSOZYME, 'rb') as fin, gzip.open(gz_file, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        obj = OpenmmSimulation(str(gz_file), {'outdir': self.outdir / 'sim'})
        obj.load_pdb()
        self.assertGreater(obj.pdb.topology.getNumAtoms(), 0)

    def test_relax_pdb(self):
        sim_outdir = self.outdir / 'relax'
        sim = OpenmmSimulation(self.fixed_pdb, {
            'outdir': sim_outdir,
            'pdb_id': 'lysozyme',
            'max_iterations': 10,
            'steps': 0,     # skip dynamics, energy minimization only
            'force_fields': ['amber14/protein.ff14SB.xml'],
        })
        info = sim.relax_pdb()
        self.assertEqual(info['status'], 'succeed')
        self.assertIn('original_energy', info)
        self.assertIn('relaxed_energy', info)
        self.assertTrue((sim_outdir / 'relaxed.pdb').exists())
        self.assertTrue((sim_outdir / 'relaxed_state.xml').exists())
        self.assertTrue((sim_outdir / 'relaxed_info.json').exists())
