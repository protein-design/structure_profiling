'''
references:
wwPDB OneDep validation reports (X-ray + EM guides, updated 2025) 
and MolProbity (Chen et al. 2010, Read et al. 2011).

Protein geometry QC for a PDB file.
1. Chain continuity
    missing backbone atoms
    chain breaks
2. Bond geometry
        Ca-Ca distance, 
        peptide C-N distance,
        backbone bond angles (N-CA-C, CA-C-N)
3. Torsion angles:
    phi/psi/omega in degrees, Ramachandran
    classification,
    cis/trans peptide bonds
4. Steric clashes:
    non-bonded heavy-atom contacts; bonded pairs
    (same / sequential residue, disulfides) and
    hydrogens are excluded, clash score is
    normalized per 1000 heavy atoms
5. Compactness           radius of gyration

execution:
    run()   -> summary statistics dict ('details' holds raw data)
    check() -> (pass/fail, reason)
    
    '''

import warnings

import numpy as np
import pandas as pd
from pathlib import Path

from ..structure import Pdb
from ..qc import Clashes
from .geometry_distance import GeometryDistance
from .geometry_angle import GeometryAngle

class Geometry(Pdb):
    # steric clash QC thresholds (canonical values live in Clashes)
    clash_score_cutoff:float = Clashes.clash_score_cutoff
    severe_clash_distance:float = Clashes.severe_clash_distance
    # MolProbity clashscore bands (vdw_clash_score): <5 very good,
    # 5-10 acceptable, 10-20 investigate, >20 poor
    vdw_clash_cutoff:float = Clashes.vdw_clash_cutoff


    def __init__(self, pdb_file:Path, is_backbone:bool=False):
        super().__init__(pdb_file)
        # is_backbone=True skips torsion-angle QC for backbone-only models;
        self.is_backbone = is_backbone
   
    def bond_length_angle_z(self, chain) -> pd.DataFrame:
        '''
        Combined covalent-geometry Z table: bond rows from
        GeometryDistance.bond_length_z plus angle rows from
        GeometryAngle.bond_angle_z, labelled by kind ('bond' / 'angle').

        One row per measured parameter: res_pos, resname, kind, atoms,
        observed, ideal, sd, z, outlier (|Z| > 5, wwPDB criterion;
        Z_OUTLIER_CUT=5.0 in both GeometryDistance and GeometryAngle).
        '''
        bond = GeometryDistance(chain).bond_length_z()
        angle = GeometryAngle(chain).bond_angle_z()
        cols = ['res_pos', 'resname', 'kind', 'atoms', 'observed',
                'ideal', 'sd', 'z', 'outlier']
        frames = []
        if not bond.empty:
            bond = bond.copy()
            bond.insert(2, 'kind', 'bond')
            frames.append(bond)
        if not angle.empty:
            angle = angle.copy()
            angle.insert(2, 'kind', 'angle')
            frames.append(angle)
        if not frames:
            return pd.DataFrame(columns=cols)
        return pd.concat(frames, ignore_index=True)[cols]

    def chain_rmsz(self, chain) -> dict:
        '''
        RMS-Z of backbone bond lengths and bond angles for one chain
        (wwPDB criterion: RMS-Z in [0, 1] is expected; |Z| > 5 are outliers).

        Cutoffs: RMS-Z itself is a magnitude (good models ~0.5-1.0);
        the former BACKBONE_GATES used bond_rmsz=2.0 / angle_rmsz=1.5
        as design tolerances; outlier counting uses Z_OUTLIER_CUT=5.0.

        Returns {'bond_rmsz', 'angle_rmsz', 'n_bond_outliers',
        'n_angle_outliers'}; RMS-Z is None when the chain has no measurable
        parameters.
        '''
        df = self.bond_length_angle_z(chain)
        if df.empty:
            return {
                'bond_rmsz': None,
                'angle_rmsz': None,
                'n_bond_outliers': 0,
                'n_angle_outliers': 0
            }
        rmsz = lambda z: float(np.sqrt(np.mean(np.square(z)))) if len(z) else None
        return {
            'bond_rmsz': rmsz(df.loc[df['kind'] == 'bond', 'z'].to_numpy()),
            'angle_rmsz': rmsz(df.loc[df['kind'] == 'angle', 'z'].to_numpy()),
            'n_bond_outliers': int(((df['kind'] == 'bond') & df['outlier']).sum()),
            'n_angle_outliers': int(((df['kind'] == 'angle') & df['outlier']).sum()),
        }

    def run(self, chain_no:str) -> dict:
        '''compute all geometry metrics of one chain.

        Thresholds applied downstream (see check()): missing>0 fails,
        ca max>5 A fails, peptide max>5 A (mean>3 A) fails, torsion null
        fraction>0.05 fails, clash_score>clash_score_cutoff (30.0) fails.

        Returns a dict of summary statistics; raw per-residue data
        (arrays / DataFrames) is nested under 'details'.
        '''
        chain = self.get_chain(chain_no)
        if chain is None:
            return {}
        geo_stat, geo_data = {}, {}

        # distance metrics: Ca-Ca, peptide C-N, radius of gyration
        gd = GeometryDistance(chain)
        geo_stat['total_residues'] = len(gd.standard_residues())
        
        # missing backbone atoms: N, CA, C, O; missing residues are skipped
        missing = gd.missing_residues()
        geo_stat['num_missing'] = len(missing)
        geo_data['missing_residues'] = missing

        ca_dist = gd.ca_distances()
        if ca_dist is not None and len(ca_dist) > 0:
            geo_stat['ca_distance'] = pd.Series(ca_dist).describe().to_dict()
            geo_data['ca_distance'] = ca_dist

        pep_dist = gd.peptide_distances()
        if pep_dist is not None and len(pep_dist) > 0:
            geo_stat['peptide_distance'] = pd.Series(pep_dist).describe().to_dict()
            geo_data['peptide_distance'] = pep_dist

        # chain continuity: counts in the summary (num_breaks / num_fragments
        # / num_strained, also the Breaks.GATES keys),
        # per-pair localisation under 'details'
        cont = gd.chain_continuity()
        geo_stat['chain_continuity'] = {k: cont[k] for k in ('num_residues', \
            'num_fragments', 'num_breaks', 'num_strained')}
        geo_data['chain_continuity'] = cont

        rg = gd.radius_of_gyration()
        if rg is not None:
            geo_stat['radius_of_gyration'] = round(rg, 2)

        clashes = gd.steric_clashes()
        if clashes is not None:
            heavy = len(gd.heavy_atoms())
            severe = clashes[clashes['distance'] < gd.severe_clash_distance] \
                if not clashes.empty else clashes
            geo_stat['num_clashes'] = len(clashes)
            geo_stat['num_severe_clashes'] = len(severe)
            geo_stat['clash_score'] = \
                round(1000.0 * len(clashes) / heavy, 2) if heavy else 0.0
            if not clashes.empty:
                geo_stat['clashes'] = clashes['distance'].describe().to_dict()
                geo_data['clashes'] = clashes
                        
        # angle metrics: bond angles, torsion angles, omega (cis/trans)
        ga = GeometryAngle(chain)
        
        bond_ang = ga.bond_angles()
        if bond_ang is not None and len(bond_ang) > 0:
            geo_stat['bond_angles'] = {
                col: pd.Series(bond_ang[col].dropna()).describe().to_dict()
                for col in ('angle_N_CA_C', 'angle_CA_C_N')
            }
            geo_data['bond_angles'] = bond_ang

        angles = ga.torsion_angles()
        if angles is not None and len(angles) > 0:
            perc = angles['ramachandran'].value_counts(normalize=True)
            geo_stat['torsion_angles'] = perc.round(2).to_dict()
            geo_data['torsion_angles'] = angles
            # omega: cis / trans peptide bonds, cis non-proline is rare
            conf = angles['conformation'].value_counts(normalize=True)
            cis = angles[angles['conformation'] == 'cis']
            geo_stat['omega'] = {
                'cis': round(float(conf.get('cis', 0.0)), 2),
                'trans': round(float(conf.get('trans', 0.0)), 2),
                'cis_nonpro': int((cis['residue'] != 'PRO').sum()),
            }

        geo_stat['details'] = geo_data
        return geo_stat

    def check(self, geo_info:dict) -> tuple:
        '''gate on summary statistics returned by run().

        Cutoffs (fail -> reason):
          num_missing > 0 -> 'missing_residues' (no tolerance);
          chain_continuity num_breaks > 0 -> 'chain_breaks'
            (num_breaks = 0);
          ca max > 5 A -> 'ca_distance' (ideal CA-CA ~3.8 A);
          peptide max > 5 A or mean > 3 A -> 'peptide_distance'
            (ideal C-N = 1.329 A);
          torsion null fraction > 0.05 -> 'torsion_angles'
            (backbone-only models skip this gate);
          num_severe_clashes > 0 -> 'severe_clashes'
            (contacts < severe_clash_distance = 2.0 A;
            former BACKBONE_GATES num_servere_clashes = 0);
          clash_score > clash_score_cutoff (30.0/1000 atoms) -> 'clashes'.
      
        '''
        if not geo_info:
            return False, 'null_geometry'
        
        # no missing atoms in any residues
        if geo_info.get('num_missing', 0):
            return False, 'missing_residues'
        # chain continuity: any real break (broken / numbering-gap peptide)
        # splits the chain; num_fragments == num_breaks + 1 (see run()),
        # gate: num_breaks = 0
        cont = geo_info.get('chain_continuity') or {}
        if cont.get('num_breaks', 0) > 0:
            return False, 'chain_breaks'
        # good cutoff 3.8A, allowed 5A
        ca = geo_info.get('ca_distance')
        if ca is not None and pd.notna(ca.get('max')) and ca['max'] > 5:
            return False, 'ca_distance'
        # ideal C-N 1.329 A; fail on a break (max > 5 A, the Breaks
        # pep_max_dist cutoff) or stretched peptides (mean > 3 A)
        pep = geo_info.get('peptide_distance')
        if pep is not None and pd.notna(pep.get('max')):
            if pep['max'] > 5 or (pd.notna(pep.get('mean')) and pep['mean'] > 3):
                return False, 'peptide_distance'
        if self.is_backbone is False:
            # percentage of Ramachandran null entries (missing phi/psi)
            torsion = geo_info.get('torsion_angles') or {}
            if torsion.get('null', 0) > 0.05:
                return False, 'torsion_angles'
        # severe steric clashes: non-bonded heavy-atom contacts below
        # severe_clash_distance (2.0 A); num_servere_clashes = 0
        if geo_info.get('num_severe_clashes', 0) > 0:
            return False, 'severe_clashes'
        # steric clashes: clash score = clashes per 1000 heavy atoms,
        # see Geometry.clash_score_cutoff
        if geo_info.get('clash_score', 0) > self.clash_score_cutoff:
            return False, 'clashes'
        return True, 'pass'

