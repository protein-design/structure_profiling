#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 sp chain_torsions ::: {0..258}
    - sp chain_torsions 0
function:
    - backbone torsion / bond-angle metrics of every chain in chain_pdb
      (phi/psi/omega with cis-trans conformation, Ramachandran,
      bond-angle Z-scores, wwPDB peptide / main-chain planarity),
      collected into table chain_torsions (src/sql/qc.sql)

Cutoffs live in geometry/torsions.py (Torsions):
    - Z_OUTLIER_CUT = 5.0     |Z| above this is a bond-angle outlier
    - CIS_CUT = 90 deg        |omega| below this is a cis peptide
    - OMEGA_DEV_CUT = 30 deg  omega dev above this is a planarity outlier
    - MAINCHAIN_OOP_CUT = 10 deg  out-of-plane above this is an outlier
    - GATES: angle_rmsz = 1.5, z_outlier_frac = 0.01,
      torsion_null_frac = 0.05, peptide_outlier_frac = 0.05,
      mainchain_outlier_frac = 0.05, rama_outlier_frac = 0.05
'''
from src.helper import *
from structure_profiling import Pdb, Torsions, TABLE_COLUMNS

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.chain_id, c.chain_no, a.relative_pdb
        FROM chain_pdb          a
        LEFT JOIN pdb_chainid   c ON a.chain_id = c.chain_id
        LEFT JOIN pdb_chain_seq s ON a.chain_id = s.chain_id
        LEFT JOIN meta_pdb      m ON a.pdb_id = m.pdb_id
        WHERE m.pdb_group={pdb_group}
            AND a.relative_pdb IS NOT NULL
            AND s.pro_len > 3
            AND c.model_no = 0
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE a.chain_id = b.chain_id
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} query files ready...")
    meta["rows"] += len(rows)
    wrapped = []
    for row in rows:
        pdb_file = Path(params['output_pdb_dir']) / row['relative_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            wrapped.append(row)
    return wrapped

def torsion_metrics(pdb_file, chain_no):
    '''Torsions.run() summary of one chain (torsions / bond angles).

    Returns backbone bond-angle stats (n_angles, angle_rmsz,
    n_z_outliers, z_outlier_frac), phi/psi/omega stats (n_torsions,
    torsion_null_frac, frac_trans, frac_cis, n_cis_nonpro),
    MolProbity-style Ramachandran (rama_analysed, rama_favoured_frac,
    rama_allowed_frac, rama_outlier_frac) and wwPDB peptide /
    main-chain planarity (n_peptides, peptide_outlier_frac,
    mainchain_outlier_frac).
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    return Torsions(chain).run(keep_details=False)

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = torsion_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: torsion metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_torsions',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
