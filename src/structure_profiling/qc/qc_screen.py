'''
prompt:
    Add unit test cases for new input data 10BL_0_A_0.pdb and 10BL_0_B_0.pdb
    into tests/test_qc_screen.py, or corresponding test files.
    Verify bond_rmsz and angle_rmsd values.
prompt :
    explain "z_outlier_frac" in the context of the code provided.
    - All returns are in markdown format
    - return formula if possible
    - recommend a threshold.
'''
from pathlib import Path

from ..structure import Pdb
from ..geometry.geometry_distance import GeometryDistance
from ..geometry.geometry_angle import GeometryAngle
from ..geometry.geometry import Geometry

  
# illustrative gates
# tune for my denovo design pipeline
# allow design tolerance
BACKBONE_GATES = {
    # hard failures
    'num_missing': 0,   # no residue may lack N/CA/C/O' atoms tolerated
    'num_breaks': 0,        # real discontinuities, broken or numbering-gap peptides fail, 
    'num_fragments': 1,     # one contiguous chain

    
    # tolerance failures (design tolerance)
    'num_strained': 3,      # stretched-but-intact peptides tolerated (tune)
    'num_clashes': 10,  # non-bonded heavy atom - heavy atom distances < 2.4 A
    'num_servere_clashes': 0, # < 2.0 A
    'clash_score': 10,  # non-bonded contacts per 1000 heavy atoms
    'angle_rmsz': 1.5, # |Z| RMS over angles, expects [0, 1]; 
    'bond_rmsz': 2.0,   # |Z| RMS over bonds, wwPDB expects [0, 1]; 
    'ca_max_dist': 5.0, # CA(i)-CA(i+1) max (A); ideal ~3.8, as check()
    'pep_max_dist': 5.0, # C(i)-N(i+1) max (A) <-- peptide C-N break; ideal 1.329 (was 5.0)
    'z_outlier_frac': 0.01, # fraction of |Z|>5 bond/angle outliers (was fixed count n_z_outliers=10)
    'peptide_outlier_frac': 0.05, # omega dev > 30 deg (was 0.05)
    'mainchain_outlier_frac': 0.05, # >10 deg out-of-plane (was 0.05)
}

# full-atom screen gates: BACKBONE_GATES (backbone) + side chains
STRUCTURE_GATES = {
    **BACKBONE_GATES,
    'cb_outlier_frac': 0.05,    # cbeta dev > CB_DEVIATION_CUT (0.25 A)
    'n_chirality_outliers': 0,  # D / planar CA centres (L-amino acids)
    'rama_outlier_frac': 0.05,  # MolProbity-style Rama outliers
    'vdw_clash_score': 20.0,  # per 1000 heavy atoms
}


# ===========================================================================
# Cheap backbone-only QC screen (needs only N/CA/C/O)
# ===========================================================================

def backbone_qc_screen(pdb_file: Path, chain_no: str = None,
        gates: dict = None, keep_details: bool = True) -> dict:
    '''Screen many backbone models with N/CA/C/O-only metrics.
    Cutoffs (defaults = BACKBONE_GATES)
    Per file: missing backbone atoms, chain continuity (breaks /
    fragments / strained peptides), CA-CA and peptide C-N distances,
    chain_rmsz (bond + angle Z vs ideal), wwPDB planarity outliers and
    the steric-clash score.
    Returns one row per (file, chain): pass, reason and metric values;
    geometry DataFrames are nested under 'details' when keep_details.
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    info = {'file': str(pdb_file), 'chain': chain_no, 'details': None,}
    if chain is None:
        info.update({
            'pass': False,
            'reason': 'no_chain',
        })
        return info
    # initialize gates with defaults if not provided; allow overrides
    gates = dict(BACKBONE_GATES) if gates is None \
        else dict(BACKBONE_GATES, **gates)

    gd = GeometryDistance(chain)
    missing = len(gd.missing_residues())
    ca = gd.ca_distances()
    pep = gd.peptide_distances()
    ca_max_dist = float(ca.max()) if len(ca) else 0.0
    pep_max_dist = float(pep.max()) if len(pep) else 0.0

    # chain continuity: numbering gaps, broken / strained peptides
    # (num_breaks / num_fragments / num_strained match the gate names)
    cont = gd.chain_continuity()
    num_breaks = cont['num_breaks']
    num_fragments = cont['num_fragments']
    num_strained = cont['num_strained']

    ga = GeometryAngle(chain)
    pl = ga.peptide_planarity()
    n_pep = len(pl)
    pep_frac = float(pl['peptide_outlier'].sum() / n_pep) if n_pep else 0.0
    mc_frac = float(pl['mainchain_outlier'].sum() / n_pep) if n_pep else 0.0
    
    geom = Geometry(pdb_file)
    rmsz = geom.chain_rmsz(chain)
    zdf = geom.bond_length_angle_z(chain)
    n_z = int(zdf['outlier'].sum()) if not zdf.empty else 0
    z_outlier_frac = float(n_z / len(zdf)) if not zdf.empty else 0.0

    # steric clashes: same metric and gate as Geometry.run()/check()
    clashes = gd.steric_clashes()
    num_clashes = len(clashes) if clashes is not None else 0
    servere_clashes = clashes[clashes['distance'] < 2.0] if clashes is not None else None
    num_servere_clashes = len(servere_clashes) if servere_clashes is not None else 0
    heavy = len(gd.heavy_atoms())
    clash_score = round(1000.0 * num_clashes / heavy, 2) if heavy else 0.0

    failures = []
    if missing > gates['num_missing']:
        failures.append('missing_residues')
    if num_breaks > gates['num_breaks']:
        failures.append('chain_breaks')
    if num_fragments > gates['num_fragments']:
        failures.append('chain_fragments')
    if num_strained > gates['num_strained']:
        failures.append('chain_strain')
    if ca_max_dist > gates['ca_max_dist']:
        failures.append('ca_distance')
    if (rmsz['bond_rmsz'] or 0) > gates['bond_rmsz']:
        failures.append('bond_rmsz')
    if (rmsz['angle_rmsz'] or 0) > gates['angle_rmsz']:
        failures.append('angle_rmsz')
    if z_outlier_frac > gates['z_outlier_frac']:
        failures.append('z_outliers')
    if pep_frac > gates['peptide_outlier_frac']:
        failures.append('peptide_planarity')
    if mc_frac > gates['mainchain_outlier_frac']:
        failures.append('mainchain_planarity')
    if pep_max_dist > gates['pep_max_dist']:
        failures.append('peptide_distance')
    if clash_score > gates['clash_score']:
        failures.append('clashes')
    
    info.update({
        'pass': not failures,
        'reason': 'pass' if not failures else ','.join(failures),
        'num_missing': missing,
        'num_breaks': num_breaks,
        'num_fragments': num_fragments,
        'num_strained': num_strained,
        'ca_max_dist': round(ca_max_dist, 2),
        'pep_max_dist': round(pep_max_dist, 2),
        'angle_rmsz': rmsz['angle_rmsz'],
        'bond_rmsz': rmsz['bond_rmsz'],
        'n_z_outliers': n_z,
        'z_outlier_frac': round(z_outlier_frac, 3),
        'peptide_outlier_frac': round(pep_frac, 3),
        'mainchain_outlier_frac': round(mc_frac, 3),
        'num_clashes': num_clashes,
        'num_servere_clashes': num_servere_clashes,
        'clash_score': clash_score,
    })
    if keep_details:
        info['details'] = {
            'bond_angle_z': zdf,
            'planarity': pl,
            'chain_breaks': cont['breaks'],
            'chain_strained': cont['strained']
        }
    return info


# ===========================================================================
# Full-atom QC screen (backbone + side chains)
# ===========================================================================

def structure_qc_screen(pdb_file:Path, chain_no: str = None,
        gates: dict = None, keep_details: bool = True) -> dict:
    '''Screen full-atom models: backbone + side-chain geometry.

    The backbone half is delegated to backbone_qc_screen() (same
    backbone metrics, gates and reason labels, incl. chain continuity);
    side-chain / full-atom gates (defaults =
    STRUCTURE_GATES): cb_outlier_frac=0.05 (fraction of residues
    with Cb deviation > CB_DEVIATION_CUT=0.25 A; Gly has no CB and cannot
    deviate); n_chirality_outliers=0 (D / planar CA centres in an
    L-amino-acid chain; |chiral volume| < CHIRAL_VOL_MIN=0.7 A^3 is
    planar); rama_outlier_frac=0.05 (residue-type aware Ramachandran
    outliers from rama_molprobity()); vdw_clash_score=
    20.0 (MolProbity vdW overlaps >= VDW_OVERLAP_CUT=0.40 A per 1000
    heavy atoms; >20 is the poor band, hydrogens excluded).

    Hetero residues (glycan NAG, waters, ligands) are not part of the
    protein chain and are excluded from every metric; a chain holding no
    standard amino acids (e.g. a glycan-only NAG chain) reports
    'no_protein'.
    Per file: the backbone_qc_screen metrics (incl. chain continuity:
    breaks / fragments / strained peptides) plus Cb deviation, CA
    chirality, Ramachandran outliers and the vdW-overlap clash score;
    raw tables are nested under 'details' when keep_details.
    Returns one row per (file, chain): pass, reason and metric values.
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    info = {'file': str(pdb_file), 'chain': chain_no, 'details': None,}
    if chain is None:
        info.update({
            'pass': False,
            'reason': 'no_chain',
        })
        return info

    gd = GeometryDistance(chain)
    if len(gd.standard_residues()) == 0:
        # chain exists but holds no standard amino acids
        # (e.g. a glycan-only NAG chain): nothing to QC
        info.update({
            'pass': False,
            'reason': 'no_protein',
        })
        return info

    gates = dict(STRUCTURE_GATES) if gates is None \
        else dict(STRUCTURE_GATES, **gates)

    # backbone half: same metrics, gates and reason labels as
    # backbone_qc_screen() (BACKBONE_GATES is a subset of STRUCTURE_GATES)
    info = backbone_qc_screen(pdb_file, chain_no=chain_no,
        gates={k: gates[k] for k in BACKBONE_GATES},
        keep_details=keep_details)

    # side chains: Cb deviation from the idealized position
    cb = gd.cbeta_deviation()
    n_cb = len(cb)
    cb_frac = float(cb['outlier'].sum() / n_cb) if n_cb else 0.0

    # CA chirality: D / planar centres in an L-amino-acid chain
    ga = GeometryAngle(chain)
    ch = ga.ca_chirality()
    n_chiral_bad = int(ch['handedness'].isin(['D', 'planar']).sum()) \
        if not ch.empty else 0

    # residue-type aware Ramachandran (analytic boxes)
    rama = ga.rama_molprobity()

    # MolProbity vdW overlaps: radii-aware; hydrogens excluded and the
    # score normalized per 1000 heavy atoms (consistent with clash_score)
    heavy = len(gd.heavy_atoms())
    vdw = gd.vdw_clashes(include_hydrogen=False)
    vdw_score = gd.vdw_clash_score(vdw, n_atoms=heavy)['clash_score'] \
        if heavy else 0.0

    # side-chain failures are appended after the backbone ones
    failures = [] if info['pass'] else info['reason'].split(',')
    if cb_frac > gates['cb_outlier_frac']:
        failures.append('cbeta_deviation')
    if n_chiral_bad > gates['n_chirality_outliers']:
        failures.append('chirality')
    if rama['outlier'] > gates['rama_outlier_frac']:
        failures.append('rama_outliers')
    if vdw_score > gates['vdw_clash_score']:
        failures.append('vdw_clashes')

    info.update({
        'pass': not failures,
        'reason': 'pass' if not failures else ','.join(failures),
        'cb_outlier_frac': round(cb_frac, 3),
        'n_chirality_outliers': n_chiral_bad,
        'rama_outlier_frac': rama['outlier'],
        'vdw_clash_score': vdw_score,
    })
    if keep_details:
        info['details'].update({
            'cbeta_deviation': cb,
            'ca_chirality': ch,
            'vdw_clashes': vdw,
            'rama_outliers': rama['outliers'],
        })
    return info