#!/usr/bin/python
'''
define funcitons used by the below
    - combo2_relaxed.py
    - combo3_relaxed.py
    - combo4_relaxed.py
    - combo5_relaxed.py
'''
from src.helper import *
from structure_profiling import WriteTmp, ProcessJson, OpenmmSimulation

@header
def retrieve_pdb(params, meta):
    # clean null data in the table
    validate_data('relaxed_pdb', params, meta, params['simulate_dir'])

    pdb_group = params['pdb_group']
    combo = params['combo']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.combo_id, c.fixed_pdb
        FROM combo{combo}_fixed c
        LEFT JOIN meta_pdb m ON c.pdb_id=m.pdb_id
        WHERE m.pdb_group={pdb_group}
            AND c.fixed_pdb IS NOT NULL
            AND c.combo_id NOT IN (
                SELECT combo_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] += len(rows)
    return rows

def relax(row, params, meta):
    fixed_pdb_file = params['output_pdb_dir'] / row['fixed_pdb']
    pdb_id, combo_id = row['pdb_id'], row['combo_id']
    endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / combo_id
    outdir = params['simulate_dir'] / endpoint

    file_name = 'relaxed.pdb'
    relaxed_pdb = outdir / file_name
    info_json = outdir / 'relaxed_info.json'
    if params['overwrite'] or not relaxed_pdb.is_file() or not info_json.is_file():
        with WriteTmp(params['tmp_dir'], outdir, params['table_name']) as writer:
            # relax pdb
            curr_params = {
                'pdb_id': pdb_id,
                'outdir': writer.temp_dir_path,
            }
            curr_params.update(params)
            info = OpenmmSimulation(fixed_pdb_file, curr_params).relax_pdb()
            meta['relax_' + info['status']] += 1
    else:
        meta['relax_skip'] += 1                    

    # load relax info
    info = ProcessJson.load_json(info_json, {})
    # update record
    relative_path = str(endpoint / file_name) if relaxed_pdb.is_file() else None
    rec = (
        pdb_id,
        combo_id,
        relative_path,
        info.get('original_energy'),
        info.get('relaxed_energy'),
    )
    if params['verbose']:
        print(f"\n### Relax {combo_id}, relaxed pdb={relative_path}")
    meta['insertion'] += 1
    return rec