#!/usr/bin/python
'''
requirements:
    - table label_chain
example:
    - parallel -j4 sp af_msa ::: {0..258}
    - sp af_msa 0
function:
    - MSA and build *.a3m for structure prediction
'''
import asyncio
from structure_profiling import Clean, async_concurr
from .func_predict import *

def retrieve(params, meta):
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    table_name = params['table_name']
    input_size = params['max_input_size']
    query = f"""
        SELECT a.pdb_id, a.input_id, a.query_faa, a.outdir
        FROM {params_table_name} a
        WHERE a.group_no={pdb_group}
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE a.input_id = b.input_id
            )
            AND a.input_size <= {input_size}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} query files ready...")
    meta["rows"] += len(rows)
    wrapped = []
    for row in rows:
        row['a3m_name'] = f"{row['input_id']}.a3m"
        row['a3m_file'] = Path(row['outdir']) / row['a3m_name']
        wrapped.append(row)
    return wrapped

# NO asynchronous funciton is needed
def post_func(input):
    # input format is consistent with returns of 
    # run_colabfold_search and run_colabfold_msa
    row, params, meta = input
    if row.get('a3m_file') and row['a3m_file'].is_file():
        # insert
        rec = (row['input_id'], row['a3m_name'], params['host'],)
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records([tuple(rec),], params['table_name'])
        meta['records'] += 1
    else:
        meta['no_a3m'] += 1


if __name__ == "__main__":
    params.update({
        'tmp_dir': params['tmp_dir'],
        'colabfold_bin_dir': os.getenv('colabfold_bin_dir'),
        'colabfold_db': os.getenv('colabfold_db'),
        'overwrite':False,
        'params_table_name': 'af_params',
        'table_name': 'af_msa',
        'max_input_size': 800,
    })
    add_params(params)

    # TODO: need robust local MSA
    # kill the mmseqs leftovers of interrupted runs (no live
    # colabfold_search driver), keep running MSA jobs untouched
    # killed = Clean.clean_msa(verbose=params['verbose'])
    # meta['clean_pids'] = len(killed)

    validate_file_rows('a3m_name', params, meta)
    
    # asynchronous functions
    funcs = [
        run_colabfold_msa, # MSA API
        # run_colabfold_search, # local
    ]
    rows = retrieve(params, meta)
    input = [(row, params, meta) for row in rows]
    asyncio.run(
        async_concurr.main(funcs, input, post_func)
    )
    
    footer(meta, params)
