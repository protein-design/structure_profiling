#!/usr/bin/python
'''
requirements:
    - table chain_pdb
    - rawdata_dir env var (original complex PDBs, meta_pdb.relative_pdb)
example:
    - parallel -j8 sp chain_ligand ::: {0..258}
    - sp chain_ligand 0
function:
    - ligand / binding-site neighborhood of every chain in chain_pdb.
      Default ligand mode reads the ORIGINAL complex from rawdata_dir:
      the split chain PDB holds no HETATM (ChainAtomSelect drops
      ligands/waters for dssp/fix/relax). Sites are the HETATM ligands
      of this chain (waters skipped), candidates the standard residues
      of this chain (Neighbors, chain_no scoped); or an explicit
      activity site via params['site_residues'] (e.g. "A/57,A/102",
      read from the chain PDB). Every residue with a heavy atom within
      search_cutoff of a site is a neighbor - neighbors may or may not
      be continuous in the chain sequence (pocket), n_fragments counts
      the sequence-continuous segments. Summary -> table chain_ligand
      (src/sql/ligand.sql), full contact table -> neighbors.csv /
      neighbors.pkl below output_pdb_dir (relative paths recorded in
      neighbor_csv / neighbor_pkl).

Cutoffs live in structure/neighbors.py (Neighbors):
    - search_cutoff = 5.0 A   heavy-atom distance flagged as a contact
    - waters (hetflag 'W', HOH / DOD / WAT) are skipped as sites
'''
from src.helper import *
from structure_profiling import Pdb, Neighbors, ProcessPickle, TABLE_COLUMNS

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT a.pdb_id, a.chain_id, c.chain_no, c.model_no,
            a.relative_pdb, m.relative_pdb AS relative_complex
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
    # default ligand mode reads the ORIGINAL complex from rawdata_dir:
    # the split chain PDB holds no HETATM (ChainAtomSelect drops
    # ligands/waters). Activity-site mode only needs the protein and
    # stays on the chain PDB. Missing files are skipped (no row) so
    # they can be picked up on a later run instead of being recorded
    # as permanent zero rows.
    site_residues = parse_site_residues(params.get('site_residues'))
    wrapped = []
    for row in rows:
        if site_residues:
            pdb_file = Path(params['output_pdb_dir']) / row['relative_pdb']
        else:
            pdb_file = Path(params['rawdata_dir']) / row['relative_complex']
        if pdb_file.is_file():
            row['pdb_file'] = pdb_file
            wrapped.append(row)
        else:
            meta['missing_input'] += 1
    return wrapped


def ligand_neighbors(pdb_file, chain_no=None, cutoff=None,
        site_residues=None, model_no=None):
    '''Neighbors.run() summary of one PDB file (None when the file is
    missing or unparsable). model_no restricts to one model of a
    multi-model structure (e.g. NMR complexes): the other models are
    detached so the residue pool is not duplicated per model.'''
    p = Pdb(pdb_file, None, False)
    if not p.structure:
        return None
    if model_no is not None:
        for model in list(p.structure):
            if model.id != model_no:
                p.structure.detach_child(model.id)
    nb = Neighbors(p.structure, chain_no=chain_no)
    return nb.run(cutoff=cutoff, site_residues=site_residues,
        keep_details=True)


def parse_site_residues(spec):
    '''"A/57,A/102" -> [('A', 57), ('A', 102)]; None passthrough.'''
    if not spec:
        return None
    pool = []
    for item in re.split(r'[,;]+', str(spec)):
        item = item.strip()
        if not item:
            continue
        chain_no, _, resseq = item.rpartition('/')
        pool.append((chain_no, int(resseq)))
    return pool


def export_neighbors(info, outdir, params, meta):
    '''export the full neighbor contact table below outdir:
    neighbors.csv (pandas) and neighbors.pkl (summary dict).

    Skipped when the file exists and not params['overwrite']. Returns
    the (csv, pkl) file names present on disk.
    '''
    outdir.mkdir(parents=True, exist_ok=True)
    # 1. contact table -> csv
    csv_name = 'neighbors.csv'
    csv_file = outdir / csv_name
    if params['overwrite'] or not csv_file.is_file():
        info['neighbors'].to_csv(csv_file, index=False)
        meta['csv_succeed'] += 1
    else:
        meta['csv_skip'] += 1
    # 2. summary -> pickle (the csv holds the contact table)
    pkl_name = 'neighbors.pkl'
    pkl_file = outdir / pkl_name
    if params['overwrite'] or not pkl_file.is_file():
        summary = {k: v for k, v in info.items() if k != 'neighbors'}
        ProcessPickle(pkl_file, params['verbose']).save_pickle(
            {'summary': summary})
        meta['pkl_succeed'] += 1
    else:
        meta['pkl_skip'] += 1
    return (csv_name if csv_file.is_file() else None,
            pkl_name if pkl_file.is_file() else None)


def run(rows, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    site_residues = parse_site_residues(params.get('site_residues'))
    for row in rows:
        try:
            info = ligand_neighbors(row['pdb_file'], row.get('chain_no'),
                cutoff=params.get('search_cutoff'),
                site_residues=site_residues,
                model_no=row.get('model_no'))
            if info is None:
                meta['null_structure'] += 1
                continue
            # export the contact table below output_pdb_dir
            endpoint = (Path(params['table_name']) / row['pdb_id'][:2]
                / row['pdb_id'] / row['chain_id'])
            outdir = params['output_pdb_dir'] / endpoint
            if info['n_sites'] > 0:
                csv_name, pkl_name = export_neighbors(
                    info, outdir, params, meta)
            else:
                csv_name = pkl_name = None
                meta['no_site'] += 1
            info.update({
                'pdb_id': row['pdb_id'],
                'chain_id': row['chain_id'],
                'neighbor_csv': str(endpoint / csv_name)
                    if csv_name else None,
                'neighbor_pkl': str(endpoint / pkl_name)
                    if pkl_name else None,
            })
            rec = [info.get(i) for i in cols]
            if params['verbose']:
                print(rec)
            yield tuple(rec)
        except Exception as e:
            print(f"WARNING: ligand neighbors of {row['pdb_file']}"
                  f" failed. error={e}")

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'table_name': 'chain_ligand',
        'search_cutoff': None,     # default 5.0 A (Neighbors.search_cutoff)
        'site_residues': None,     # activity site, e.g. "A/57,A/102"
    })
    add_params(params)

    rows = retrieve(params, meta)
    record_iter = run(rows, params, meta)
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta, params)
