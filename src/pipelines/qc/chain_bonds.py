#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 sp chain_bonds ::: {0..258}
    - sp chain_bonds 0
function:
    - covalent bond-length metrics of every chain in chain_pdb (RMS-Z vs
      the Engh & Huber ideal bond lengths, |Z| > 5 outliers), collected
      into table chain_bonds (src/sql/qc.sql)

Cutoffs live in qc/geometry_distance.py (GeometryDistance):
    - Z_OUTLIER_CUT = 5.0   |Z| above this is a covalent-geometry outlier
    - IDEAL_BONDS           Engh & Huber (mean, sd) targets: N-CA, CA-C,
                            C-O, CA-CB plus the C(i)-N(i+1) peptide bridge
    - GATES (former qc_screen BACKBONE_GATES): bond_rmsz = 2.0 (wwPDB expects
      RMS-Z in [0, 1]), z_outlier_frac = 0.01 (bond + angle rows combined;
      this table holds the bond-length half)
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

def bond_metrics(pdb_file, chain_no):
    '''Covalent bond-length summary of one chain.

    Z-scores of the backbone bond lengths (N-CA, CA-C, C-O, CA-CB and
    the C(i)-N(i+1) peptide bridge) vs the Engh & Huber ideals in
    GeometryDistance.IDEAL_BONDS. Returns n_bonds, bond_rmsz (RMS of
    the Z-scores, the kind='bond' half of Geometry.chain_rmsz; wwPDB
    expects [0, 1], former BACKBONE_GATES tolerates 2.0), n_bond_outliers
    (|Z| > Z_OUTLIER_CUT = 5), bond_z_outlier_frac and bond_max_absz.
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    bonds = GeometryDistance(chain).bond_length_z()
    n_bonds = len(bonds)
    if not n_bonds:
        return {
            'n_bonds': 0,
            'bond_rmsz': None,
            'n_bond_outliers': 0,
            'bond_z_outlier_frac': 0.0,
            'bond_max_absz': 0.0,
        }
    z = bonds['z'].to_numpy()
    n_out = int(bonds['outlier'].sum())
    return {
        'n_bonds': n_bonds,
        'bond_rmsz': float(np.sqrt(np.mean(np.square(z)))),
        'n_bond_outliers': n_out,
        'bond_z_outlier_frac': round(n_out / n_bonds, 3),
        'bond_max_absz': round(float(np.abs(z).max()), 3),
    }

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = bond_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: bond metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_bonds',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
