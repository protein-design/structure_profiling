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
    # dssp
    'chain_dssp': ['pdb_id', 'chain_id', 'host', 'model_no',
        'chain_no', 'dssp_seq',],
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
