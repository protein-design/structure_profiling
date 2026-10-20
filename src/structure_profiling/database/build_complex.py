'''
insert data into empty table
Subset of bio_omics bioomics.database.build_complex: only the tables
written by the structure profiling pipelines are kept.
'''
import os
from .insert_db import InsertDb
from .delete_db import DeleteDb

# build SQL query string
TABLE_COLUMNS = {
    # split complex into chain combos (combo2..5)
    'combo2': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo3': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo4': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    'combo5': ['pdb_id', 'combo_id', 'chain_combo', 'chain_pdbs',
        'combo_pdb', 'combo_csv', 'combo_fasta'],
    # dssp
    'chain_dssp': ['pdb_id', 'chain_id', 'host', 'model_no',
        'chain_no', 'dssp_seq',],
    # fix
    'chain_fixed': ['pdb_id', 'chain_id', 'fixed_pdb',],
    'combo2_fixed': ['pdb_id', 'combo_id', 'fixed_pdb',],
    'combo3_fixed': ['pdb_id', 'combo_id', 'fixed_pdb',],
    'combo4_fixed': ['pdb_id', 'combo_id', 'fixed_pdb',],
    'combo5_fixed': ['pdb_id', 'combo_id', 'fixed_pdb',],
    # relax
    'chain_relaxed': ['pdb_id', 'chain_id', 'relaxed_pdb',
        'original_energy', 'relaxed_energy',],
    'combo2_relaxed': ['pdb_id', 'combo_id', 'relaxed_pdb',
        'original_energy', 'relaxed_energy',],
    'combo3_relaxed': ['pdb_id', 'combo_id', 'relaxed_pdb',
        'original_energy', 'relaxed_energy',],
    'combo4_relaxed': ['pdb_id', 'combo_id', 'relaxed_pdb',
        'original_energy', 'relaxed_energy',],
    'combo5_relaxed': ['pdb_id', 'combo_id', 'relaxed_pdb',
        'original_energy', 'relaxed_energy',],
    #qc
    'pdb_qc': ['pdb_id', 'chain_id', 'angle_rmsz', 'bond_rmsz',
        'ca_max_dist', 'cb_outlier_frac', 'clash_score',
        'mainchain_outlier_frac', 'n_chirality_outliers',
        'n_z_outliers', 'num_breaks', 'num_clashes', 'num_fragments',
        'num_missing', 'num_strained', 'num_servere_clashes',
        'pep_max_dist', 'peptide_outlier_frac', 'rama_outlier_frac',
        'vdw_clash_score', 'z_outlier_frac',],
    'pdb_relaxed_qc': ['pdb_id', 'chain_id', 'angle_rmsz', 'bond_rmsz',
        'ca_max_dist', 'cb_outlier_frac', 'clash_score',
        'mainchain_outlier_frac', 'n_chirality_outliers',
        'n_z_outliers', 'num_breaks', 'num_clashes', 'num_fragments',
        'num_missing', 'num_strained', 'num_servere_clashes',
        'pep_max_dist', 'peptide_outlier_frac', 'rama_outlier_frac',
        'vdw_clash_score', 'z_outlier_frac',],
    # clashes (pipelines/qc/chain_clashes.py, src/sql/qc.sql)
    'chain_clashes': ['pdb_id', 'chain_id', 'n_heavy_atoms',
        'num_clashes', 'num_severe_clashes', 'clash_score',
        'n_vdw_clashes', 'vdw_clash_score', 'max_overlap',],
    # chain breaks / continuity (pipelines/qc/chain_breaks.py, src/sql/qc.sql)
    'chain_breaks': ['pdb_id', 'chain_id', 'n_residues',
        'num_missing', 'num_breaks', 'num_fragments', 'num_strained',
        'ca_max_dist', 'pep_max_dist',],
    # meta
    'meta': ['pdb_group', 'host', 'ip_addr', 'script'],
    'meta_execute': ['table_name', 'pdb_group', 'host', 'ip_addr', 'script'],
}


class BuildComplex(InsertDb, DeleteDb):

    def __init__(self, verbose:bool=False, chunk_size:int=None):
        '''
        args: data is iterable object of pdb_file, or record
        '''
        InsertDb.__init__(self, verbose, chunk_size)
        DeleteDb.__init__(self, verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def insert_batch_records(self, data, table_name:str, cols=None):
        if cols is None:
            cols = TABLE_COLUMNS[table_name]
        fmt = ['%s',] * len(cols)
        cols_str = ', '.join(cols)
        fmt_str= ', '.join(fmt)
        insert_query = f"""
            INSERT INTO {table_name} ({cols_str})
            VALUES ({fmt_str})
        """
        # load insertion
        func = lambda rec: [rec,]
        if self.verbose:
            print(f"\nBegin to insert records into {table_name}...", end='')
        return self.func_insert(data, insert_query, func)
