'''
functions:
    - used by af_rmsd.py, af_tmalign.py
Decoupled from bio_omics: bioomics imports are replaced by structure_profiling.
'''

from src.helper import *
from structure_profiling import EvalPred

def retrieve_pred(params, meta):
    params_table_name = params['params_table_name']
    table_name = params['pred_table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.input_id, a.pred_outdir
        FROM {table_name}             a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        WHERE p.group_no = {pdb_group}
            AND a.pred_outdir IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['predictions'] = len(rows)
    print(f"Detect {len(rows)} predictions in {table_name}...")
    return rows

def pred_summary(params):
    params_table = params['params_table_name']
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.relaxed_pdb_name AS pdb_name
        FROM {table_name}        a
        LEFT JOIN {params_table} p ON a.input_id = p.input_id
        LEFT JOIN meta_pdb       m ON p.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND a.relaxed_pdb_name IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    relaxed_pdb_names = [i['pdb_name'] for i in rows]
    print(f"Detect {len(rows)} predicted pdb files in {table_name}...")
    return relaxed_pdb_names


def chain_raw_model_pdb(params):
    pdb_status = 'raw'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['output_pdb_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.relative_pdb AS native_pdb
        FROM {summary_table_name}     a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN chain_pdb           c ON a.input_id = c.chain_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.relative_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.relative_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, 'chain', query

def chain_fixed_model_pdb(params):
    pdb_status = 'fixed'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['output_pdb_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.fixed_pdb AS native_pdb
        FROM {summary_table_name}     a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN chain_fixed         c ON a.input_id = c.chain_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.fixed_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.fixed_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, 'chain', query

def chain_relaxed_model_pdb(params):
    pdb_status = 'relaxed'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['simulate_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.relaxed_pdb AS native_pdb
        FROM {summary_table_name}     a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN chain_relaxed       c ON a.input_id = c.chain_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.relaxed_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.relaxed_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, 'chain', query

def combo_raw_model_pdb(params, combo):
    pdb_status = 'raw'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['output_pdb_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.combo_pdb AS native_pdb
        FROM {summary_table_name}        a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN {combo}             c ON a.input_id = c.combo_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.combo_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.combo_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, combo, query

def combo_fixed_model_pdb(params, combo):
    pdb_status = 'fixed'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['output_pdb_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.fixed_pdb AS native_pdb
        FROM {summary_table_name}     a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN {combo}_fixed       c ON a.input_id = c.combo_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.fixed_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.fixed_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, combo, query
 
def combo_relaxed_model_pdb(params, combo):
    pdb_status = 'relaxed'
    pdb_group = params['pdb_group']
    params_table_name = params['params_table_name']
    summary_table_name = params['summary_table_name']
    table_name = params['table_name']
    native_outdir = params['simulate_dir']
    query = f"""
        SELECT p.input_id, '{pdb_status}' AS pdb_status, p.pdb_id, p.category, 
            p.outdir AS pred_outdir, a.relaxed_pdb_name AS pred_pdb_name,
            '{native_outdir}' AS native_outdir, c.relaxed_pdb AS native_pdb
        FROM {summary_table_name}     a
        LEFT JOIN {params_table_name} p ON a.input_id = p.input_id
        LEFT JOIN {combo}_relaxed     c ON a.input_id = c.combo_id
        WHERE p.group_no = {pdb_group}
            AND a.ranking = 1
            AND c.relaxed_pdb IS NOT NULL
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE c.relaxed_pdb = b.native_pdb
            )
    ;"""
    return pdb_status, combo, query

def calculate_rmsd(rows:list, meta):
    input_rows = []
    for row in rows:
        row['pdb1'] = Path(row['native_outdir']) / row['native_pdb']
        row['pdb2'] = Path(row['pred_outdir']) / row['pred_pdb_name']
        if row['pdb1'].is_file() and row['pdb2'].is_file():
            input_rows.append(row)
    meta['pred_pdb_files'] = len(input_rows)
    
    rmsd_row_iter = EvalPred().batch_rmsd(input_rows)
    for row in rmsd_row_iter:
        rec = (
            row['input_id'],
            row['pdb_status'],
            row['rmsd'],
            row['native_pdb'],
            row['pred_pdb_name'],
        )
        print(rec)
        yield rec
        
def calculate_tmalign(rows:list, params, meta):
    prefix = f"{params['pdb_status']}_{params['category']}"
    ep = EvalPred(params['mmalign_dir'])
    
    for row in rows:
        native_pdb = Path(row['native_outdir']) / row['native_pdb']
        model_pdb = Path(row['pred_outdir']) / row['pred_pdb_name']
        if native_pdb.is_file() and model_pdb.is_file():
            try:
                res = ep.TMalign(native_pdb, model_pdb)
                if res:
                    meta[f'{prefix}_succeed'] += 1
                else:
                    meta[f'{[prefix]}_fail'] += 1
                rec = (
                    row['input_id'],
                    row['pdb_status'],
                    res.get('RMSD'),
                    res.get('TM1'),
                    res.get('TM2'),
                    row['native_pdb'],
                    row['pred_pdb_name'],
                )
                print(rec)
                yield rec
            except Exception as e:
                meta[f'{prefix}_fail'] += 1
                pass