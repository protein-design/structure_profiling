#!/usr/bin/python
'''
requirements:
    - table chain_fixed
example:
    - parallel -j8 sp chain_fixed_breaks ::: {0..258}
    - sp chain_fixed_breaks 0
function:
    - chain-continuity metrics of every chain in chain_fixed (missing
      backbone atoms, breaks, fragments, strained peptides), collected
      into table chain_fixed_breaks (src/sql/qc.sql)

Cutoffs live in geometry/breaks.py (Breaks):
    - cn_break = 2.0 A       C-N above this is a real break (ideal 1.329 A)
    - ca_break = 4.5 A       CA-CA above this is a break (ideal ~3.8 A)
    - cn_strain = 1.5 A      C-N above this but <= cn_break is strained
    - GATES: num_missing = 0, num_breaks = 0, num_fragments = 1,
      num_strained = 3, ca_max_dist = 5.0, pep_max_dist = 5.0
'''
from src.helper import *
from structure_profiling import Pdb, Breaks, TABLE_COLUMNS

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.chain_id, c.chain_no, a.fixed_pdb
        FROM chain_fixed        a
        LEFT JOIN pdb_chainid   c ON a.chain_id = c.chain_id
        LEFT JOIN pdb_chain_seq s ON a.chain_id = s.chain_id
        LEFT JOIN meta_pdb      m ON a.pdb_id = m.pdb_id
        WHERE m.pdb_group={pdb_group}
            AND a.fixed_pdb IS NOT NULL
            AND s.pro_len > 3
            AND c.model_no = 0
            AND NOT EXISTS (
                SELECT 1 FROM {table_name} b
                WHERE a.chain_id = b.chain_id
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} query files ready...")
    meta["rows"] += len(rows)
    wrapped = []
    for row in rows:
        pdb_file = Path(params['output_pdb_dir']) / row['fixed_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            wrapped.append(row)
    return wrapped

def break_metrics(pdb_file, chain_no):
    '''Breaks.run() summary of one chain (continuity / breaks).

    Returns n_residues, num_missing (residues lacking N/CA/C/O),
    num_breaks (broken or numbering-gap peptides), num_fragments
    (= 1 + num_breaks), num_strained (1.5 A < C-N <= 2.0 A) plus
    ca_max_dist and pep_max_dist.
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    return Breaks(chain).run(keep_details=False)

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = break_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: break metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_fixed_breaks',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
