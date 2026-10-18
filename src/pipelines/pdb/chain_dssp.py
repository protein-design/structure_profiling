'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j8 bomics chain_dssp ::: {0..241}
    - bomics chain_dssp 0
functions:
    - Analyze secondary structure by chains
    - put dssp sequence to table chain_dssp
'''

from src.helper import *
from structure_profiling import AnalyzeDssp, TABLE_COLUMNS

@header
def search_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.chain_id, c.chain_pdb
        FROM chain_pdb          c
        LEFT JOIN meta_pdb      m ON c.pdb_id = m.pdb_id
        LEFT JOIN pdb_chain_seq s ON c.chain_id = s.chain_id
        LEFT JOIN pdb_chainid   i ON c.chain_id = i.chain_id
        WHERE m.pdb_group = {pdb_group}
            AND s.chain_seq IS NOT NULL
            AND i.model_no = 0
            AND c.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = params['output_pdb_dir'] / row['chain_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            yield row

def pdb_dssp(row_iter, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in row_iter:
        item = {
            'host': params['host'],
            'pdb_id': row['pdb_id'],
            'chain_id': row['chain_id'],
        }
        try:
            # note: pdb_file contains only one chain
            ad = AnalyzeDssp(row['pdb_file'])
            for model, dssp_data in ad.dssp_iter():
                item['model_no'] = model.id
                for chain_no, dssp in dssp_data.items():
                    item['chain_no'] = chain_no
                    item['dssp_seq'] = ''.join([i['dssp_code'] for i in dssp])
                    meta['records'] += 1
                    break
                break
        except Exception as e:
            print(f"ERROR: Fail in pdb_secondary(). error={e}")
            meta['fail'] += 1
        rec = [item.get(col) for col in cols]
        yield tuple(rec)


if __name__ == "__main__":
    params.update({
        'chunk_size': 10,
        'table_name': 'chain_dssp',
    })
    row_iter = search_pdb(params, meta)
    record_iter = pdb_dssp(row_iter, params, meta)

    print(f"Try to analyze secondary structure and insert results into {params['table_name']}...")
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)