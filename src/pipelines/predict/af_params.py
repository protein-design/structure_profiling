#!/usr/bin/python
'''
requirements:
    - table label_chain
example:
    - parallel -j8 sp af_params ::: {0..258}
    - sp af_params 0
function:
    - build params for pdb structure prediction
'''
from src.helper import *
from structure_profiling import ProcessFasta
from .func_predict import *

def build(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        query_faa = params['input_dir'] / row['relative_faa']
        if query_faa.is_file():
            input_seq = ProcessFasta(query_faa).get_first_seq()
            meta[f"{params['combo']}_records"] += 1
            outdir = params['outdir'] / pdb_id[:2] /pdb_id / row['input_id']
            rec = (
                pdb_id,
                row['input_id'],
                params['input_table_name'],
                params['input_col_name'],
                len(input_seq),
                params['pdb_group'],
                params['combo'],
                row['categories'],
                str(query_faa),
                str(outdir),
            )
            yield rec
        else:
            meta[f"{params['combo']}_skip"] += 1

if __name__ == "__main__":
    params.update({
        'table_name': 'af_params',
    })
    add_params(params)

    # AF2
    qc = QueryComplex(params['verbose'])
    pool = (
        retrieve_chain(params, 'chain_af'),
        retrieve_combo(params, 'combo2_af', '2'),
        retrieve_combo(params, 'combo3_af', '3'),
        retrieve_combo(params, 'combo4_af', '4'),
        retrieve_combo(params, 'combo5_af', '5'),
    )
    for curr_params, query in pool:
        rows = qc.list_data(query)
        print(f"{len(rows)} PDB files ready...")
        meta[f"{curr_params['outdir_name']}_rows"] += len(rows)        

        record_iter = build(rows, curr_params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records(record_iter, params['table_name'])
    
    footer(meta, params)
