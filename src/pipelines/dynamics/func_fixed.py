#!/usr/bin/python
'''
define funcitons used by the below
    - combo2_fixed.py
    - combo3_fixed.py
    - combo4_fixed.py
    - combo5_fixed.py
'''
from src.helper import *
from structure_profiling import ProcessJson, OpenmmPdbfixer

@header
def retrieve_pdb(params, meta):
    # clean null data in the table
    validate_data('fixed_pdb', params, meta)

    q = QueryComplex(params['verbose'])
    rows = q.retrieve_combo(params['pdb_group'], params['combo'], params['table_name'])
    print(f"pdb_group={params['pdb_group']}, pdb_files={len(rows)}")
    meta['rows'] += len(rows)
    return rows

def fix(rows, params, meta):
    for row in rows:
        pdb_id, combo_id = row['pdb_id'], row['combo_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / combo_id
        outdir = params['output_pdb_dir'] / endpoint

        # fix pdb
        file_name = 'fixed.pdb'
        fixed_pdb = outdir / file_name
        info_json = outdir / 'fixed_info.json'
        if params['overwrite'] or not fixed_pdb.is_file() or not info_json.is_file():
            input_pdb = params['output_pdb_dir'] / row['combo_pdb']
            curr_params = {
                'verbose': params['verbose'],
                'outdir': outdir,
            }
            op = OpenmmPdbfixer(input_pdb, curr_params)
            info = op.fast_fixer()
            meta['fix_' + info['status']] += 1
        else:
            meta['fix_skip'] += 1

        # load fix_info
        info = ProcessJson.load_json(info_json, {})

        # update record
        relative_path = str(endpoint / file_name) if fixed_pdb.is_file() else None
        rec = (pdb_id, combo_id, relative_path)
        if params['verbose']:
            print('Record: ', rec)
        meta['insertion'] += 1
        yield rec
