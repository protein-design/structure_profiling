'''
MD simulation using openMM
Trajactory file and relaxed pdb could be used for further analysis,
such as RMSD, docking simulation
'''
from datetime import datetime
import gzip
import openmm as mm
from openmm import app
import openmm.unit as unit
import os
from pathlib import Path
from sys import stdout

from .timeout import timeout
from ..common import ProcessJson

FORCE_FIELDS = {
    # a comprehensive set of Amber14 force field files
    'general': ['amber14-all.xml',],
    'protein': ['amber14/protein.ff14SB.xml', 'amber99sb.xml',],
    'membrane_protein': ['charmm36.xml',],
    'nucleic': ['amber14/all_nucleic.xml',],
    'DNA': ['amber14/DNA.OL15.xml',],
    'RNA': ['amber14/RNA.OL3.xml',],
    # explicit solvent: water
    'explicit_solvent': [
        'amber14/tip3p.xml',    # most common, standard TIP3P
        'amber14/tip3pfb.xml',  # updated version of tip3p
        'amber14/spce.xml',     # SPC/E
        'amber14/opc.xml',      # higher accuracy
    ],
    'implicit_solvent': ['implicit/gbn2.xml', 'implicit/obc.xml',],
    'small_molecule': ['gaff.xml', 'gaff2.xml', 'charmm36/chem.xml',],
    # combined with explicit solvent
    'ion': ['amber14/ions.xml',],
}

# water system
COMPLEX_FF = {
    'general_use': ['amber14-all.xml', 'amber14/tip3p.xml',],
    'protein_protein': ['amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml',],
    'membrane_protein': ['charmm36.xml', 'amber14/tip3pfb.xml',],
    'channel_protein': ['charmm36.xml', 'amber14/ions.xml', 'amber14/tip3pfb.xml',],
    'protein_dna': ['amber14/protein.ff14SB.xml', 'amber14/DNA.OL15.xml'],
    'protein_rna': ['amber14/protein.ff14SB.xml', 'amber14/RNA.OL13.xml'],
    'protein_ion': ['amber14/protein.ff14SB.xml', 'amber14/ions.xml', 'amber14/tip3pfb.xml',],
    'protein_ligand': ['amber14/protein.ff14SB.xml', 'gaff.xml', 'amber14/tip3pfb.xml',],
}

class OpenmmSimulation:

    def __init__(self, pdb_file:str, args:dict):
        '''
        note: openmm only accept string path rather than Path object
        '''
        self.pdb_file = str(pdb_file)
        self.args = args
        self.outdir = str(self.args.get('outdir'))
        Path(self.outdir).mkdir(parents=True, exist_ok=True)
        self.verbose = self.args.get('verbose', False)
        self.pdb = None
        self.simulation = None 
        # load_pdb() records the periodic-box provenance here, so the
        # dict must exist even when load_pdb() is called on its own
        self.info = {
            'pdb_file': self.pdb_file,
            'pdb_id': self.args.get('pdb_id'),
        }

    @timeout(30)
    def relax_pdb(self):
        '''
        note: max duration < 10 min in larges-scale screening
        '''
        self.info = {'pdb_id': self.args.get('pdb_id'),}
        if self.verbose:
            print(f"\nTry to relax {self.pdb_file}:")
        start_time = datetime.now()

        try:
            # simulation pipeline
            self.load_pdb()
            self.create_system()
            self.minimize_energy()
            self.export_state()
            self.info['status'] = 'succeed'
            
            # calculate relaxed energy
            minmized_state = self.simulation.context.getState(getEnergy=True)
            self.info["relaxed_energy"] = str(minmized_state.getPotentialEnergy())

            delta = datetime.now() - start_time
            self.info['duration_min'] = round(delta.total_seconds() / 60)
            if self.verbose:
                print(f"Relaxing {self.pdb_file} is done | {datetime.now()}")
        except TimeoutError as e:
            self.info['timeout_relax_pdb'] = str(e)
            self.info['status'] = 'timeout'
        except Exception as e:
            self.info['fail_relax_pdb'] = str(e)
            self.info['status'] = 'fail'

        # save self.info
        if self.verbose:
            print('### Relax info:', self.info)
        if os.path.isdir(self.outdir):
            json_file = os.path.join(self.outdir, 'relaxed_info.json')
            ProcessJson.save_json(self.info, json_file)
        return self.info

    def load_pdb(self):
        '''
        create self.pdb
        '''
        if self.pdb_file.endswith('.gz'):
            with gzip.open(self.pdb_file, 'rt') as f:
                # app.PDBFile(file, ...) has no 'pdbfile' keyword
                self.pdb = app.PDBFile(f)
        else:
            self.pdb = app.PDBFile(self.pdb_file)
        self.info['periodic_box'] = 'input_pdb'
        
        # Add box dimensions (if they don't exist or are a dummy):
        # OpenMM's PDBFile puts a dummy 0.1 nm box on PDB files without a
        # CRYST1 record; PME with a 1 nm cutoff needs at least a 2 nm box,
        # else create_system() raises "cutoff distance cannot be greater
        # than half the periodic box size"
        box = self.pdb.topology.getPeriodicBoxVectors()
        if box is None or (
                box[0][0] + box[1][1] + box[2][2]) < 1.5 * unit.nanometer:
            if self.verbose:
                print(f"Add periodic box vector for data from {self.pdb_file}")
            self.pdb.topology.setPeriodicBoxVectors([
                mm.Vec3(2.0, 0.0, 0.0),
                mm.Vec3(0.0, 2.0, 0.0),
                mm.Vec3(0.0, 0.0, 2.0),
            ] * unit.nanometer)
            self.info['periodic_box'] = 'default'

    def create_system(self):
        '''
        create self.simulation
        '''
        # define force field
        ff = self.args.get('force_fields', ['amber99sb.xml', 'tip3p.xml'])
        self.info['force_field'] = ff
        forcefield = app.ForceField(*ff)

        # initialize molecular system
        system = forcefield.createSystem(
            self.pdb.topology,
            nonbondedMethod = app.PME, # PME periodic, Nocutoff
            nonbondedCutoff = 1 * unit.nanometer,
            constraints = app.HBonds, # HBonds, AllBonds, None
            rigidWater=True,        # Keep water molecules rigid
            hydrogenMass = 1.5 * unit.amu,
        )

        # Langevin dynamics, with the LFMiddle discretization
        # (J. Phys. Chem. A 2019, 123, 28, 6056-6079).
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, # 300K is room temperature
            1.0 / unit.picosecond, # friction coefficient
            0.004 * unit.picoseconds, # 1-4 fs for stability
            # 2.0 * unit.femtoseconds,
        )

        # simulator
        try:
            cuda_platform = mm.Platform.getPlatformByName('CUDA')
            platformProperties = {"Precision": "single"}
            self.simulation = app.Simulation(
                self.pdb.topology, system, integrator,
                cuda_platform, platformProperties
            )
            if self.verbose:
                print(f"Try to minimize energy of {self.pdb_file} using CUDA...")
        except Exception as e:
            print(f"Warning: openmm can't load CUDA, error={e}.")
            self.simulation = app.Simulation(
                self.pdb.topology, system, integrator
            )
            if self.verbose:
                print(f"Try to minimize energy of {self.pdb_file} using CPU...")
        self.simulation.context.setPositions(self.pdb.positions)
        # original energy
        curr_state = self.simulation.context.getState(getEnergy=True)
        self.info['original_energy'] = str(curr_state.getPotentialEnergy())

    def minimize_energy(self):
        '''
        update self.simulation
        '''
        # energy minimization
        max_iterations = self.args.get('max_iterations', 0)
        self.simulation.minimizeEnergy(maxIterations=max_iterations)

        # define exports
        steps = self.args.get('steps', 10_000)
        report_interval = 1_000
        if os.path.isdir(self.outdir):
            # checkpoint
            outfile = os.path.join(self.outdir, 'checkpoint.chk')
            self.simulation.reporters.append(
                app.CheckpointReporter(outfile, report_interval)
            )
            # dynamic coordinate data
            outfile = os.path.join(self.outdir, 'trajactory.dcd')
            self.simulation.reporters.append(
                app.DCDReporter(outfile, report_interval)
            )
            # pdb at each state
            outfile = os.path.join(self.outdir, 'structure.pdb')
            self.simulation.reporters.append(
                app.PDBReporter(outfile, steps)
            )
            # meta data of simulation
            outfile = os.path.join(self.outdir, 'state_data.csv')
            self.simulation.reporters.append(
                app.StateDataReporter(outfile, report_interval, step=True, \
                    potentialEnergy=True, temperature=True)
            )
        else:
            self.simulation.reporters.append(
                app.StateDataReporter(stdout, report_interval, step=True, \
                    potentialEnergy=True, temperature=True)
            )
        # simulation
        self.simulation.step(steps)

    def export_state(self):
        '''
        export self.simulation
        '''
        # Save the current state of the simulation to a file.
        outfile = os.path.join(self.outdir, 'relaxed_state.xml')
        self.simulation.saveState(outfile)

        # Save minimized structure
        relaxed_pdb_file = os.path.join(self.outdir, 'relaxed.pdb')
        self.info['relaxed_pdb'] = relaxed_pdb_file
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        with open(relaxed_pdb_file, 'w') as f:
            app.PDBFile.writeFile(self.simulation.topology, positions, f)