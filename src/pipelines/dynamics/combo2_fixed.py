#!/usr/bin/python
'''
requirements:
    - table combo2
example:
    - parallel -j8 sp combo2_fixed ::: {0..258}
    - sp combo2_fixed 0
function:
    - fix pdb files
    - insert data into table combo2_fixed
'''
from src.helper import *
from .func_fixed import *


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'combo': 2,
        'table_name': 'combo2_fixed',
    })
    
    rows = retrieve_pdb(params, meta)
    record_iter = fix(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)