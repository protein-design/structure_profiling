#!/usr/bin/python
'''
requirements:
    - tables af_msa
example: 
    - parallel -j4 sp af_pred ::: {0..241}
    - sp af_pred 0
functions:
    - predict structure
'''

from src.helper import *
from .func_predict import *

if __name__ == "__main__":
    params.update({
        'tmp_dir': params['tmp_dir'],
        'colabfold_bin_dir': os.getenv('colabfold_bin_dir'),
        'relax': True,
        'params_table_name': 'af_params',
        'msa_table_name': 'af_msa',
        'table_name': 'af_pred',
        'max_input_size': 800,
    })
    
    # a3m files
    rows = retrieve_msa(params, meta)

    bc = BuildComplex(params['verbose'])
    for row in rows:
        rec = run_colabfold_batch(row, params, meta)
        if rec:
            bc.insert_batch_records([rec,], params['table_name'])

    footer(meta)
