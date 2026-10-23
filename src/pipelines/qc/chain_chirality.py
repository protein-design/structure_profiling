#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 sp chain_chirality ::: {0..258}
    - sp chain_chirality 0
function:
    - CA-chirality metrics of every chain in chain_pdb (handedness of
      every CA centre via the signed N-CA-C-CB chiral volume), collected
      into table chain_chirality (src/sql/qc.sql)

Cutoffs live in qc/geometry_angle.py (GeometryAngle):
    - CHIRAL_VOL_L = 2.6 A^3    ideal L-amino-acid CA chiral volume
    - CHIRAL_VOL_MIN = 0.7 A^3  0 < volume < this is a 'planar' centre;
      volume < 0 is a D-amino acid (rare); glycin (no CB) is 'NA'
    - GATE (former qc_screen STRUCTURE_GATES): n_chirality_outliers = 0
      (D / planar CA centres in an L-amino-acid chain)
'''
from src.helper import *
from structure_profiling import Pdb, GeometryAngle, TABLE_COLUMNS

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

def chirality_metrics(pdb_file, chain_no):
    '''CA-chirality summary of one chain (stereo validation).

    Returns n_residues (standard residues), n_chirality_analysed
    (residues with a CB; glycin is 'NA'), n_chirality_outliers
    ('D' + 'planar' centres; former STRUCTURE_GATES gate 0) and
    chirality_outlier_frac (outliers / analysed residues).
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    ch = GeometryAngle(chain).ca_chirality()
    n_res = len(ch)
    analysed = ch[ch['handedness'] != 'NA'] if n_res else ch
    n_out = int(analysed['handedness'].isin(['D', 'planar']).sum()) \
        if len(analysed) else 0
    return {
        'n_residues': n_res,
        'n_chirality_analysed': int(len(analysed)),
        'n_chirality_outliers': n_out,
        'chirality_outlier_frac': round(n_out / len(analysed), 3)
            if len(analysed) else 0.0,
    }

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = chirality_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: chirality metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_chirality',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
