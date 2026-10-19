#!/usr/bin/python
'''
requirements:
    - The table label_complex is ready.
example:
    - parallel -j8 sp combo5 ::: {0..258}
    - sp combo5 0
function:
    - split PDB of complex to paired chains PDB.
    - export complex csv and fasta for structure predication
note:
    - only consider model.id == 0
'''
from src.helper import *
from .func_combo import *


##################################
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'combo': 5,
        'table_name': 'combo5',
    })

    row_iter = retrieve_data(params, meta)
    record_iter = build_data(row_iter, params, meta)

    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)