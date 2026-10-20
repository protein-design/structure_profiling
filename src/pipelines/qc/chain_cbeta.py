#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 sp chain_cbeta ::: {0..258}
    - sp chain_cbeta 0
function:
    - MolProbity-style Cbeta-deviation metrics of every chain in
      chain_pdb (distance of the observed CB from the idealized
      position built on the N/CA/C frame), collected into table
      chain_cbeta (src/sql/qc.sql)

Cutoffs live in geometry/geometry_distance.py (GeometryDistance):
    - CB_BOND_LEN = 1.521 A     idealized CA-CB bond
    - CB_ANGLE = 110.0 deg      idealized N-CA-CB = C-CA-CB
    - CB_DEVIATION_CUT = 0.25 A deviation flagged as a CB outlier
    - GATE (former qc_screen STRUCTURE_GATES): cb_outlier_frac = 0.05;
      glycin (no CB) cannot deviate and reports cb_dev None
'''
from src.helper import *
from structure_profiling import Pdb, GeometryDistance, TABLE_COLUMNS

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

def cbeta_metrics(pdb_file, chain_no):
    '''Cbeta-deviation summary of one chain (side-chain placement).

    Returns n_residues (standard residues), n_cb_analysed (residues
    with a measurable CB; glycin reports cb_dev None), cb_max_dev
    (largest deviation, A) and cb_outlier_frac (residues with
    cb_dev > CB_DEVIATION_CUT = 0.25 A over standard residues, the
    same fraction as the former structure_qc_screen; gate 0.05).
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    cb = GeometryDistance(chain).cbeta_deviation()
    n_res = len(cb)
    n_out = int(cb['outlier'].sum())
    devs = cb['cb_dev'].dropna()
    return {
        'n_residues': n_res,
        'n_cb_analysed': int(len(devs)),
        'cb_max_dev': round(float(devs.max()), 3) if len(devs) else 0.0,
        'cb_outlier_frac': round(n_out / n_res, 3) if n_res else 0.0,
    }

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = cbeta_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: cbeta metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_cbeta',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
