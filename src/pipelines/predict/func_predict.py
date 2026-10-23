'''
fuctions for structure prediction
Ported from bio_omics src/pipelines/predict/func_predict.py.
'''
import asyncio
from src.helper import *
from structure_profiling import Predictor, WriteTmp


def validate_file_rows(col_name, params, meta):
    params_table = params['params_table_name']
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT a.input_id, b.outdir, a.{col_name}
        FROM {table_name}        a
        LEFT JOIN {params_table} b ON a.input_id = b.input_id
        LEFT JOIN meta_pdb       m ON b.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    del_rows = []
    for row in rows:
        path = Path(row['outdir']) / row[col_name]
        if not path.is_file():
            del_rows.append(row['input_id'])
            
    # delete null rows
    if del_rows:
        print(f"Delete {len(del_rows)} *.asm...")
        dc = DeleteComplex(params['verbose'])
        dc.remove_null_rows(table_name, 'input_id', del_rows)
        meta['delete_rows'] = len(del_rows)


def retrieve_chain(params, outdir_name:str):
    input_table_name = 'label_chain'
    input_col_name = 'chain_id'
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS(
            SELECT a.first_chain_id AS input_id, GROUP_CONCAT(DISTINCT a.category
                ORDER BY a.category SEPARATOR ',' ) AS categories
            FROM {input_table_name} a
            LEFT JOIN meta_pdb      m ON a.pdb_id = m.pdb_id
            WHERE m.pdb_group = {pdb_group}
            GROUP BY a.first_chain_id
        )
        SELECT t.input_id, t.categories, f.pdb_id, f.relative_faa
        FROM tmp t, pdb_chain_seq f
        WHERE t.input_id = f.chain_id
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} c
                WHERE c.input_id = t.input_id
            )
    """
    cat_params = dict(params)
    cat_params.update({
        'input_table_name': input_table_name,
        'input_col_name': input_col_name,
        'input_dir': params['output_pdb_dir'],
        'combo': 1,
        'outdir_name': outdir_name,
        'outdir': params['predict_dir'] / outdir_name,
    })
    return cat_params, query

def retrieve_combo(params, outdir_name:str, combo:int):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    input_table_name = 'label_complex'
    input_col_name = 'combo_id'
    query = f"""
        WITH tmp AS(
            SELECT a.combo_id AS input_id, GROUP_CONCAT(DISTINCT a.category
                ORDER BY a.category SEPARATOR ',' ) AS categories
            FROM {input_table_name} a
            LEFT JOIN meta_pdb      p ON a.pdb_id = p.pdb_id
            WHERE p.pdb_group = {pdb_group}
                AND a.combo = {combo}
                AND NOT EXISTS (
                    SELECT 1 FROM {table_name} b
                    WHERE a.combo_id = b.input_id
                )
            GROUP BY a.combo_id
        )
        SELECT t.input_id, t.categories, c.pdb_id, c.combo_fasta AS relative_faa
        FROM tmp t, combo{combo}  c
        WHERE t.input_id = c.combo_id
    """
    cat_params = dict(params)
    cat_params.update({
        'input_table_name': input_table_name,
        'input_col_name': input_col_name,
        'input_dir': params['output_pdb_dir'],
        'outdir_name': outdir_name,
        'combo': combo,
        'outdir_name': outdir_name,
        'outdir': params['predict_dir'] /outdir_name,
    })
    return cat_params, query

def build_curr_params(row, params):
    # update output_predict_dir for exporting outputs
    pdb_id, input_id = row['pdb_id'], row['input_id']
    endpoint = Path(params['outdir_name']) / pdb_id[:2] / pdb_id / input_id

    curr_params = dict(params)
    curr_params.update({
        'predicted_at': datetime.now(),
        'host': params['host'],
        'pdb_id': pdb_id,
        'input_id': input_id,
        'endpoint': endpoint,
        'faa_file': params['input_dir'] / row['relative_faa'],
        'outdir': params['predict_dir'] / endpoint,
    })
    if params['verbose']:
        print('Arguments for prediction: ', curr_params)
    return curr_params

async def run_colabfold_search(input):
    row, params, meta = input
    # check if MSA is available
    if params['overwrite'] or not row['a3m_file'].is_file():
        outdir = row['outdir']
        tmp_name = f"colabfold_search_{params['pdb_group']}_{row['input_id']}"
        with WriteTmp(params['tmp_dir'], outdir, tmp_name) as writer:
            tmp_row = dict(row)
            tmp_row['outdir'] = writer.temp_dir_path
            # run MSA
            res = await asyncio.to_thread(
                Predictor.colabfold_search,
                tmp_row,
                params,
            )
            # move the results from the ramdisk to the local disk only if
            # the MSA succeeded and the *.a3m was delivered into the
            # temporary dir, otherwise discard the failed results
            delivered = writer.temp_dir_path / row['a3m_name']
            if res == 0 and delivered.is_file():
                meta['msa_good'] += 1
            else:
                meta['msa_fail'] += 1
                writer.move_temp = False
                print(f"Skip moving failed results of {row['input_id']}: MSA failed or no *.a3m delivered.")
    else:
        meta['no_a3m'] += 1
    return row, params, meta

async def run_colabfold_msa(input):
    row, params, meta = input
    # check if MSA is available
    if params['overwrite'] or not row['a3m_file'].is_file():
        outdir = row['outdir']
        tmp_name = f"colabfold_msa_{params['pdb_group']}_{row['input_id']}"
        with WriteTmp(params['tmp_dir'], outdir, tmp_name) as writer:
            tmp_row = dict(row)
            tmp_row['outdir'] = writer.temp_dir_path
            # run MSA 
            res = await asyncio.to_thread(
                Predictor.colabfold_batch_msa,
                tmp_row,
                params
            )
            # move the results from the ramdisk to the local disk only if
            # the MSA succeeded and the *.a3m was delivered into the
            # temporary dir, otherwise discard the failed results
            delivered = writer.temp_dir_path / row['a3m_name']
            if res == 0 and delivered.is_file():
                meta['msa_good'] += 1
            else:
                meta['msa_fail'] += 1
                writer.move_temp = False
                print(f"Skip moving failed results of {row['input_id']}: MSA failed or no *.a3m delivered.")
    
            # delete or compress for saving space
            for path in tmp_row['outdir'].iterdir():
                if path.name.endswith('_env') or path.name.endswith('_pairgreedy'):
                    shutil.rmtree(path)
                elif path.name in ('cite.bibtex', 'config.json'):
                    path.unlink()
                elif path.suffix in ('.json',):
                    cmd = ['gzip', str(path),]
                    Predictor.run_subprocess(cmd)

    else:
        meta['no_a3m'] += 1
    return row, params, meta

def retrieve_msa(params, meta):
    params_table = params['params_table_name']
    msa_table = params['msa_table_name']
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    input_size = params['max_input_size']
    query = f"""
        SELECT a.input_id, a.a3m_name, b.pdb_id, b.outdir
        FROM {msa_table}         a
        LEFT JOIn {params_table} b  ON b.input_id = a.input_id
        WHERE b.group_no = {pdb_group}
            AND a.a3m_name IS NOT NULL
            AND NOT EXISTS (
                SELECT * FROM {table_name} c
                WHERE c.input_id = a.input_id
            )
            AND b.input_size <= {input_size}
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['msa_rows'] = len(rows)
    
    frows = []
    for row in rows:
        row['outdir'] = Path(row['outdir'])
        row['input_file'] = row['outdir'] / row['a3m_name']
        if row['input_file'].is_file():
            frows.append(row)
    print(f"Detect {len(frows)} *.asm, group={pdb_group}...")

    return frows


def run_colabfold_batch(row, params, meta):
    # check if prediction is done
    done_name = f"{row['input_id']}.done.txt"
    done_file = row['outdir'] / done_name
    if not done_file.is_file():
        outdir = row['outdir']
        tmp_name = f"colabfold_batch_{params['pdb_group']}_{row['input_id']}"
        with WriteTmp(params['tmp_dir'], outdir, tmp_name) as writer:
            #keys: tools_dir, input_file, outdir
            tmp_row = dict(row)
            tmp_row['outdir'] = writer.temp_dir_path
            res =  Predictor.colabfold_batch(tmp_row, params)
            if res == 0:
                meta['succeed'] += 1
            else:
                meta['fail'] += 1

            # delete or comparess for saving space
            for path in tmp_row['outdir'].iterdir():
                name = path.name
                if name in ('cite.bibtex', 'config.json'):
                    path.unlink()
                elif path.suffix in ('.json',):
                    cmd = ['gzip', str(path),]
                    Predictor.run_subprocess(cmd)
                elif name.endswith('_env') or name.endswith('_all') \
                    or name.endswith('_greedy'):
                    shutil.rmtree(path)
    else:
        meta['skip'] += 1
    
    if done_file.is_file():
        rec = (row['input_id'], str(row['outdir']), params['host'])
        print(f"Try to insert {rec} into table {params['table_name']}")
        return rec