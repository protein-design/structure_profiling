#!/usr/bin/python
'''
requirements:
    - chain is labelled in table label_chain
example:
    - parallel -j8 sp chain_fixed ::: {0..241}
    - sp chain_fixed 0
function:
    - fix pdb files
    - insert data into table chain_fixed
'''
from src.helper import *
from structure_profiling import OpenmmPdbfixer

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT C.pdb_id, C.chain_id, C.relative_pdb
        FROM chain_pdb          C
        LEFT JOIN pdb_chainid   I ON C.chain_id = I.chain_id
        LEFT JOIN pdb_chain_seq S ON C.chain_id = S.chain_id
        LEFT JOIN meta_pdb      P ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND I.model_no = 0
            AND LENGTH(S.chain_seq) > 3
            AND P.relative_pdb IS NOT NULL
            AND C.chain_id IN (
                SELECT chain_id from label_chain
                EXCEPT
                SELECT chain_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

def fix(rows, params, meta):
    for row in rows:
        pdb_id, chain_id = row['pdb_id'], row['chain_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / chain_id
        outdir = params['output_pdb_dir'] / endpoint

        # fix pdb
        file_name = 'fixed.pdb'
        fixed_pdb = outdir / file_name
        if params['overwrite'] or not fixed_pdb.is_file():
            input_pdb = params['output_pdb_dir'] / row['relative_pdb']
            curr_params = {
                'verbose': params['verbose'],
                'outdir': outdir,
            }
            op = OpenmmPdbfixer(input_pdb, curr_params)
            fix_info = op.fast_fixer()
            meta['fix_' + fix_info['status']] += 1
        else:
            meta['fix_skip'] += 1

        # update record
        relative_path = str(endpoint / file_name) if fixed_pdb.is_file() else None
        rec = (pdb_id, chain_id, relative_path)
        if params['verbose']:
            print('Record: ', rec)
        yield rec


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'chain_fixed',
    })
    print("Try to null data in the table...")
    validate_data('fixed_pdb', params, meta)
    
    rows = retrieve_pdb(params, meta)
    record_iter = fix(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
