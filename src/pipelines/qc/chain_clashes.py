#!/usr/bin/python
'''
requirements:
    - table chain_pdb
example:
    - parallel -j8 bomics chain_clashes ::: {0..258}
    - bomics chain_clashes 0
function:
    - steric / MolProbity vdW clash metrics of every chain in chain_pdb,
      collected into table chain_clashes (src/sql/qc.sql)

Cutoffs live in structure/clashes.py (Clashes):
    - clash_cutoff = 2.4 A          non-bonded heavy-atom contacts below
    - severe_clash_distance = 2.0 A contacts below are severe clashes
    - clash_score_cutoff = 30.0     max steric clashes / 1000 heavy atoms
    - vdw_clash_cutoff = 20.0       MolProbity clashscore poor band
    - VDW_OVERLAP_CUT = 0.40 A      vdW overlap flagged as a clash
'''
from src.helper import *
from structure_profiling import Pdb, Clashes, TABLE_COLUMNS

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.chain_id, c.chain_no, a.relative_pdb
        FROM chain_pdb          a
        LEFT JOIN pdb_chainid   c ON a.chain_id = c.chain_id
        LEFT JOIN pdb_chain_seq s ON a.chain_id = s.chain_id
        LEFT JOIN meta_pdb      m ON a.pdb_id = m.pdb_id
        WHERE m.pdb_group={pdb_group}
            AND a.relative_pdb IS NOT NULL
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
        pdb_file = Path(params['output_pdb_dir']) / row['relative_pdb']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            wrapped.append(row)
    return wrapped

def clash_metrics(pdb_file, chain_no):
    '''Clashes.run() summary of one chain (steric + MolProbity vdW).

    Returns n_heavy_atoms, num_clashes (< 2.4 A), num_severe_clashes
    (< 2.0 A), clash_score (per 1000 heavy atoms), n_vdw_clashes and
    vdw_clash_score (overlaps >= 0.40 A per 1000 heavy atoms, hydrogens
    excluded) plus max_overlap.
    '''
    chain = Pdb(pdb_file).get_chain(chain_no)
    return Clashes(chain).run(keep_details=False)

def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    for row in rows:
        try:
            info = clash_metrics(row['pdb_file'], row['chain_no'])
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
                'host': params['host'],
            })
            rec = [info.get(i) for i in cols]
            print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: clash metrics of {row['pdb_file']} failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite':False,
        'table_name': 'chain_clashes',
    })
    add_params(params)
    
    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
