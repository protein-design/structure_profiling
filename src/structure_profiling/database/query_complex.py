'''
query the COMPLEX2 database
Subset of bio_omics bioomics.database.query_complex: only the queries
used by the structure profiling pipelines are kept.
'''
import os
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
