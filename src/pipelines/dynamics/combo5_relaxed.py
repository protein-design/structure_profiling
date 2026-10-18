#!/usr/bin/python
'''
requirements:
    - table combo5_fixed
example:
    - parallel -j8 bomics combo5_relaxed ::: {0..241}
    - bomics combo5_relaxed 0
function:
    - relax pdb
    - insert data into table combo5_relaxed
'''
from src.helper import *
from .func_relaxed import *


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'tmp_dir': os.getenv('tmp_dir'),
        'combo': 5,
        'table_name': 'combo5_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 2_000,
        'steps': 10_000,
    })
    
    rows = retrieve_pdb(params, meta)
    for row in rows:
        print(row)
        rec = relax(row, params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records([rec,], params['table_name'])

    footer(meta)