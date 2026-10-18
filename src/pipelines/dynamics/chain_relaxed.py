#!/usr/bin/python
'''
requirements:
    - chain is labelled in table label_chain
example:
    - parallel -j8 bomics chain_relaxed ::: {0..241}
    - bomics chain_relaxed 0
function:
    - relax pdb
    - insert data into table chain_relaxed
'''
from src.helper import *
from structure_profiling import OpenmmSimulation, ProcessJson

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT f.pdb_id, f.chain_id, f.fixed_pdb
        FROM chain_fixed   f 
        LEFT JOIN meta_pdb p ON f.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND f.fixed_pdb IS NOT NULL
            AND f.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    for row in rows:
        row['fixed_pdb_file'] = params['output_pdb_dir'] / row['fixed_pdb']
        if row['fixed_pdb_file'].is_file():
            yield row

def relax(row, params, meta):
    time0 = datetime.now()
    pdb_id, chain_id = row['pdb_id'], row['chain_id']
    endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / chain_id
    outdir = params['simulate_dir'] / endpoint

    # relax pdb
    file_name = 'relaxed.pdb'
    relaxed_pdb = outdir / file_name
    if params['overwrite'] or not relaxed_pdb.is_file():
        curr_params = {
            'pdb_id': pdb_id,
            'outdir': outdir,
        }
        curr_params.update(params)
        c = OpenmmSimulation(row['fixed_pdb_file'], curr_params)
        relax_info = c.relax_pdb()
        meta['relax_' + relax_info['status']] += 1
    else:
        meta['relax_skip'] += 1

    # load info
    info_json = outdir / 'relaxed_info.json'
    info = ProcessJson.load_json(info_json, {})

    # build record
    relative_path = str(endpoint / file_name) if relaxed_pdb.is_file() else None
    rec = (
        pdb_id,
        chain_id,
        time0,
        relative_path,
        info.get('original_energy'),
        info.get('relaxed_energy'),
    )
    if params['verbose']:
        print('Output of relaxation: ', relative_path)
    return rec


if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'chain_relaxed',
        'force_fields': ['amber14-all.xml', 'amber14/tip3pfb.xml'],
        'max_iterations': 2_000,
        'steps': 10_000,
    })
    print("Try to null data in the table...")
    validate_data('relaxed_pdb', params, meta, params['simulate_dir'])
    
    rows = retrieve_pdb(params, meta)
    for row in rows:
        rec = relax(row, params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records([rec,], params['table_name'])

    footer(meta)