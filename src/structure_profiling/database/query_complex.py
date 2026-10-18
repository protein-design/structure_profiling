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
