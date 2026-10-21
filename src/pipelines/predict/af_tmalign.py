'''
example: 
    - parallel -j8 sp af_tmalign ::: {0..241}
    - sp af_tmalign 0
functions:
    - calculate TM-align score between model and sample structure
'''

from src.helper import *
from .func_summary import *


if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'mmalign_dir': Path(os.getenv('mmalign_dir')),
        'params_table_name': 'af_params',
        'summary_table_name': 'af_summary',
        'table_name': 'af_tmalign',
    })

    pool = (
        chain_raw_model_pdb(params),
        combo_raw_model_pdb(params, 'combo2'),
        combo_raw_model_pdb(params, 'combo3'),
        combo_raw_model_pdb(params, 'combo4'),
        combo_raw_model_pdb(params, 'combo5'),
        
        chain_fixed_model_pdb(params),
        combo_fixed_model_pdb(params, 'combo2'),
        combo_fixed_model_pdb(params, 'combo3'),
        combo_fixed_model_pdb(params, 'combo4'),
        combo_fixed_model_pdb(params, 'combo5'),
        
        chain_relaxed_model_pdb(params),
        combo_relaxed_model_pdb(params, 'combo2'),
        combo_relaxed_model_pdb(params, 'combo3'),
        combo_relaxed_model_pdb(params, 'combo4'),
        combo_relaxed_model_pdb(params, 'combo5'),
    )
    for pdb_status, category, query in pool:
        params.update({
            'pdb_status': pdb_status,
            'category': category,
        })
        rows = QueryComplex(params['verbose']).list_data(query)
        prefix = f"{pdb_status}_{category}"
        meta[f'{prefix}_rows'] = len(rows)

        record_iter = calculate_tmalign(rows, params, meta)
        bc = BuildComplex(params['verbose'], params['chunk_size'])
        bc.insert_batch_records(record_iter, params['table_name'])
    footer(meta)