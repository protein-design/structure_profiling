#!/usr/bin/python
'''
define functions used by the below
    - combo2.py
    - combo3.py
    - combo4.py
    - combo5.py
note:
    - only consider model.id == 0
    - decoupled from bio_omics src/pipelines/pdb/func_combo.py
'''
from src.helper import *
from structure_profiling import Pdb, ProcessSeq


@header
def retrieve_data(params, meta):
    # clean null data in the table
    for col_name in ('combo_pdb', 'combo_csv', 'combo_fasta'):
        validate_data(col_name, params, meta)

    # retrieve data
    pdb_group = params['pdb_group']
    qc = QueryComplex(params['verbose'])
    #combo_ids
    chains_num = params['combo']
    table_name = params['table_name']
    query = f"""
        SELECT c.pdb_id, c.model_no, c.combo_id, c.chain_combo,
            c.chain_label, c.chain_ids, p.relative_pdb
        FROM label_complex c
        LEFT JOIN meta_pdb p on c.pdb_id = p.pdb_id
        WHERE p.pdb_group = {pdb_group}
            AND combo = {chains_num}
            AND NOT EXISTS (
                SELECT 1 from {table_name} a
                WHERE a.combo_id = c.combo_id
            )
    ;"""
    rows = qc.list_data(query)
    ids = [i['combo_id'] for i in rows]
    print('Check duplicates:', len(rows), len(set(ids)))
    meta[f'combo_{chains_num}'] = len(rows)

    # chain_ids
    pdb_chains = qc.pdb_chains(pdb_group)
    # add key: input_pdb, chain_seqs, chain_ids
    for row in rows:
        # absolute path
        row['input_pdb'] = params['rawdata_dir'] / row['relative_pdb']
        # parse files and sequences
        key = (row['pdb_id'], row['model_no'])
        if key in pdb_chains:
            chain_ids = row['chain_ids'].split(',')
            sub = pdb_chains[key].sort_values('chain_id')
            sub = sub[sub['chain_id'].isin(chain_ids)]
            # list type
            row['chain_seqs'] = sub['chain_seq'].to_list()
            row['chain_pdbs'] = sub['relative_pdb'].str.cat(sep=',')
            row['chain_nos'] = row['chain_combo'].split('-')
        else:
            print(f"###Warning: no match {row}")
        yield row


def build_data(row_iter, params, meta):
    '''
    build *.csv, *.fasta, *.pdb
    '''
    for row in row_iter:
        pdb_id = row['pdb_id']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
        outdir = params['output_pdb_dir'] / endpoint
        outdir.mkdir(parents=True, exist_ok=True)
        combo_id = row['combo_id']

        # 1. export combo pdb
        pdb_file_name = f"{combo_id}.pdb"
        pdb_file = outdir / pdb_file_name
        if params['overwrite'] or not pdb_file.is_file():
            p = Pdb(row['input_pdb'], None, params['verbose'])
            p.load_structure()
            status = p.to_chain_pdb(row['model_no'], row['chain_nos'], pdb_file)
            meta['combo_pdb_succeed'] += 1
        else:
            meta['combo_pdb_skip'] += 1

        # 2. export complex sequence in csv
        csv_file_name = f"{combo_id}.csv"
        csv_file = outdir / csv_file_name
        if params['overwrite'] or not csv_file.is_file():
            combo_seq = {combo_id: row['chain_seqs'],}
            status = ProcessSeq.to_colabfold_csv(combo_seq, csv_file)
            meta['combo_csv_'+ status] += 1
        else:
            meta['combo_csv_skip'] += 1

        # 3. export complex sequence in fasta
        fasta_file_name = f"{combo_id}.fasta"
        fasta_file = outdir / fasta_file_name
        if params['overwrite'] or not fasta_file.is_file():
            combo_seq = {combo_id: row['chain_seqs'],}
            status = ProcessSeq.to_colabfold_fasta(combo_seq, fasta_file)
            meta['combo_fasta_'+ status] += 1
        else:
            meta['combo_fasta_skip'] += 1

        # prepare record
        bounded_pdb = str(endpoint / pdb_file_name) if pdb_file.is_file() else None
        complex_csv = str(endpoint / csv_file_name) if csv_file.is_file() else None
        complex_fasta = str(endpoint / fasta_file_name) if fasta_file.is_file() else None
        rec = (
            pdb_id,
            combo_id,
            row['chain_combo'],
            row['chain_pdbs'],
            bounded_pdb,
            complex_csv,
            complex_fasta
        )
        if params['verbose']:
            print(params['pdb_group'], bounded_pdb)
        yield rec