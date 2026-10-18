#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 bomics pdb_qc ::: {0..241}
    - bomics pdb_qc 0
function:
    - geometry of structure
'''
from src.helper import *
from structure_profiling import structure_qc_screen, TABLE_COLUMNS

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

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = structure_qc_screen(
                row['pdb_file'], row['chain_no'], None, False
            )
            info.update(row)
            rec = [info.get(i) for i in cols]
            print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: QC of {row['pdb_file']} failed. error={e}")


if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'pdb_qc',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)