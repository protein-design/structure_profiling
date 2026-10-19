'''
query the COMPLEX2 database
Subset of bio_omics bioomics.database.query_complex: only the queries
used by the structure profiling pipelines are kept.
'''
import os

import pandas as pd

from .query_db import QueryDb

class QueryComplex(QueryDb):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def pop_pdb_group(self, table_name):
        query = f"""
            select pdb_group from meta_pdb
            except
            select pdb_group from meta_execute
            where table_name = '{table_name}'
            order by pdb_group
        ;"""
        res = self.list_data(query)
        if len(res) > 0:
            return res[0]['pdb_group']
        return None

    def retrieve_combo(self, pdb_group:int, combo:int, dest_table_name:str):
        query = f"""
            SELECT c.pdb_id, c.combo_id, c.combo_pdb
            FROM combo{combo}  c
            LEFT JOIN meta_pdb p ON c.pdb_id = p.pdb_id
            WHERE p.pdb_group = {pdb_group}
                AND c.combo_pdb IS NOT NULL
                AND c.combo_id NOT IN (
                    SELECT combo_id FROM {dest_table_name}
                )
        ;"""
        rows = self.list_data(query)
        if self.verbose:
            print(f"{len(rows)} PDB files ready...")
        return rows

    def pdb_chains(self, pdb_group:int) -> dict:
        '''
        used by the combo2..5 split pipelines
        '''
        query = f"""
            SELECT c.pdb_id, c.model_no, c.chain_id, d.relative_pdb,
                s.relative_faa AS chain_faa, s.chain_seq
            FROM pdb_chainid        c
            LEFT JOIN meta_pdb      p ON c.pdb_id=p.pdb_id
            LEFT JOIN pdb_chain_seq s ON c.chain_id=s.chain_id
            LEFT JOIN chain_pdb     d ON c.chain_id=d.chain_id
            WHERE p.pdb_group = {pdb_group}
                AND c.model_no = 0
                AND s.pro_len >= 3
                AND s.relative_faa IS NOT NULL
                AND d.relative_pdb IS NOT NULL
        ;"""
        rows = self.list_data(query)
        df = pd.DataFrame(rows)
        if not df.empty:
            g = df.groupby(['pdb_id', 'model_no'])
            return {key: sub for key, sub in g}
        return {}
