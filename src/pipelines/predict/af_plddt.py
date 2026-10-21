'''
requirements:
    - table af_summary
example: 
    - parallel -j8 sp af_plddt ::: {0..241}
    - sp af_plddt 0
functions:
    - calculate pLDDT by chain given prediction pdb
'''
from src.helper import *
from structure_profiling import ProcessPickle, CalculatePlddt

def validate_data(params):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.input_id, p.outdir, a.pred_pdb_name, a.relative_pkl
        FROM {table_name}         a
        LEFT JOIn af_params       p ON a.input_id = p.input_id
        LEFT JOIN meta_pdb        m ON p.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    del_pool = set()
    for row in rows:
        pdb_file = Path(row['outdir']) / row['pred_pdb_name']
        if not pdb_file.is_file():
            del_pool.add(row['input_id'])
        pkl_file = params['predict_dir'] / row['relative_pkl']
        if not pkl_file.is_file():
            del_pool.add(row['input_id'])
            
    #
    dc = DeleteComplex(params['verbose'])
    dc.remove_null_rows(table_name, 'input_id', list(del_pool))
    
# @header
def retrieve_rows(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT p.pdb_id, p.outdir, a.input_id, a.relaxed_pdb_name AS pred_pdb_name
        FROM af_summary           a
        LEFT JOIN af_params       p ON a.input_id = p.input_id
        LEFT JOIN meta_pdb        m ON p.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND a.relaxed_pdb_name IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} c
                WHERE a.relaxed_pdb_name = c.pred_pdb_name
            )
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['raws'] = len(rows)
    print(f"Detect {len(rows)} ready for pLDDT...")
    
    for row in rows:
        pdb_file = Path(row['outdir']) / row['pred_pdb_name']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            row['pkl_file_name'] = pdb_file.stem + '.pkl'
            meta['relaxed_pdb'] = len(rows)
            yield row

def cal_plddt(rows, params, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id / row['input_id']
        outdir = params['predict_dir'] / endpoint

        pfile = outdir / row['pkl_file_name']
        if params['overwrite'] or not pfile.is_file():
            outdir.mkdir(parents=True, exist_ok=True)
            # load structure before calculation
            cp = CalculatePlddt(row['pdb_file'], outdir)
            cp.load_structure()
            if params['verbose']:
                print(f"Try to calclate pLDDT by chain with the input {row['pdb_file']}...")
            data = cp.complex_ca_plddt()
            if data.get('pLDDT'):
                if params['verbose']:
                    print(f"Save pLDDT results to {pfile}")
                status = ProcessPickle(pfile).save_pickle(data)
                meta[status] += 1
            else:
                meta['null'] += 1
        else:
            meta['skip'] += 1

        # build record
        data = ProcessPickle(pfile).load_pickle()
        if data and data.get('pLDDT'):
            meta['insert_chains'] += 1
            for item in data['pLDDT']:
                meta['insert_records'] += 1
                rec = (
                    row['input_id'],
                    row['pred_pdb_name'],
                    str(endpoint / row['pkl_file_name']),
                    item['model'].id,
                    item['chain'].id,
                    item['stats'].get('mean_plddt'),
                )
                yield rec

 
if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'af_plddt',
    })
    validate_data(params)
        
    # retrieve pdb_ids
    rows = retrieve_rows(params, meta)

    # calculate pLDDT
    plddt_data_iter = cal_plddt(rows, params, meta)

    # insert data
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(plddt_data_iter, params['table_name'])

    footer(meta)
