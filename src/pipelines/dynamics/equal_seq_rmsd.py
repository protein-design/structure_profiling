#!/usr/bin/python
'''
requirements:
    - chain PDB files of table chain_pdb (built by bio_omics)
    - tables pdb_chain_seq, label_chain, pdb_chainid and meta_pdb
example:
    - parallel -j8 sp equal_seq_rmsd ::: {0..258}
    - sp equal_seq_rmsd 0
functions:
    - calculate CA-RMSD of chain pairs of which sequences are identical
    - put data to table equal_seq_rmsd
Decoupled from bio_omics (src/pipelines/pdb/equal_seq_rmsd.py):
the RMSD comes from the structure_profiling EvalPred (Biopython
Superimposer) instead of pyrosetta.
'''
from src.helper import *
from structure_profiling import EvalPred

def validate(params):
    table_name = params['table_name']
    pdb_group = params['pdb_group']
    query = f"""
        SELECT e.chain_id, e.chain_pdb, e.equal_chain_pdb
        FROM {table_name}     e
        LEFT JOIN pdb_chainid c ON e.chain_id = c.chain_id
        LEFT JOIN meta_pdb    m ON c.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND (e.chain_pdb IS NOT NULL OR e.equal_chain_pdb IS NOT NULL)
    """
    pool = set()
    rows = QueryComplex(params['verbose']).list_data(query)
    for row in rows:
        if row['chain_pdb']:
            file = params['output_pdb_dir'] / row['chain_pdb']
            if not file.is_file():
                pool.add(row['chain_id'])
        if row['equal_chain_pdb']:
            file = params['output_pdb_dir'] / row['equal_chain_pdb']
            if not file.is_file():
                pool.add(row['chain_id'])
    # delete records
    values = list(pool)
    if len(values):
        print(f"Try to delete invalid files: {len(values)}")
        dc = DeleteComplex(params['verbose'])
        dc.delete_rows(table_name, 'chain_id', values)

@header
def pull_chain_pdb(params, meta):
    query = """
        SELECT chain_id, relative_pdb
        FROM chain_pdb
        WHERE relative_pdb IS NOT NULL
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['num_chains'] = len(rows)
    df = pd.DataFrame(rows)
    res = dict(zip(df['chain_id'], df['relative_pdb']))
    return res

@header
def pull_equal_rmsd(params, meta):
    table_name = params['table_name']
    query = f"""
        SELECT chain_id, equal_chain_id
        FROM {table_name}
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    pairs = [set((i['chain_id'], i['equal_chain_id'])) for i in rows]
    meta['equal_pairs'] += len(rows)
    return pairs

def pull_equal_ids(pairs, params, meta):
    pdb_group = params['pdb_group']
    query = f"""
        WITH tmp AS (
            SELECT first_chain_id
            FROM label_chain   c
            LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
            WHERE p.pdb_group = {pdb_group}
        )
        SELECT chain_id, first_chain_id
        FROM pdb_chain_seq
        WHERE chain_seq IS NOT NULL
            AND first_chain_id IN (
                SELECT first_chain_id FROM tmp
            )
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    if rows:
        df = pd.DataFrame(rows)
        for _, sub in df.groupby('first_chain_id'):
            if len(sub) > 1:
                meta['num_equal_seq'] += 1
                ids = sub['chain_id'].tolist()
                combo_iter = itertools.combinations(ids, 2)
                # exclude records already in table equal_seq_rmsd
                for (a, b) in combo_iter:
                    if set((a, b)) not in pairs:
                        yield (a, b)

def calculate_rmsd(ids_pool, chain_pdb, params, meta):
    # prepare input pdbs
    pool = []
    for (id1, id2) in ids_pool:
        # chain_id is the key
        relative1 = chain_pdb.get(id1)
        relative2 = chain_pdb.get(id2)
        if relative1 and relative2:
            # pdb files
            pdb1 = params['output_pdb_dir'] / relative1
            pdb2 = params['output_pdb_dir'] / relative2
            if pdb1.is_file() and pdb2.is_file():
                rec = {
                    'id1': id1,
                    'relative1': relative1,
                    'pdb1': pdb1,
                    'id2': id2,
                    'relative2': relative2,
                    'pdb2': pdb2,
                }
                pool.append(rec)

    # calculate CA-RMSD after optimal superposition,
    # batch_rmsd yields the pairs of which rmsd is computed
    good = list(EvalPred().batch_rmsd(pool))
    meta['succeed'] = len(good)
    meta['fail'] = len(pool) - len(good)
    for rec in good:
        yield (rec['id1'], rec['id2'], rec['rmsd'],
            rec['relative1'], rec['relative2'])

if __name__ == "__main__":
    params.update({
        'chunk_size': 20,
        'table_name': 'equal_seq_rmsd',
    })
    # remove records if pdb files are missing
    validate(params)

    # prepare data
    chain_pdb = pull_chain_pdb(params, meta)
    pairs = pull_equal_rmsd(params, meta)
    ids_pool = pull_equal_ids(pairs, params, meta)

    # calculate rmsd of equal-sequence chain pairs
    record_iter = calculate_rmsd(ids_pool, chain_pdb, params, meta)

    # put data to table
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)