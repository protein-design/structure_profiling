'''
delete records of the COMPLEX2 database
Subset of bio_omics bioomics.database.delete_complex: only the deletes
used by the structure profiling pipelines are kept.
'''
import os
from .delete_db import DeleteDb

class DeleteComplex(DeleteDb):

    def __init__(self, verbose:bool=None):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def remove_null_rows(self, table_name:str, col_name:str, del_pool:list):
        n = len(del_pool)
        print(f"Try to delete {n} rows from {table_name}...")
        num_del = 0
        if n > 0:
            print(f"Number of {n} rows are invalid.", end=' ')
            chunk_size = 100
            for i in range(0, n, chunk_size):
                values = del_pool[i:i+chunk_size]
                res = self.delete_rows(table_name, col_name, values)
                if res:
                    num_del += len(values)
                    print(i, end=' ')
            print(f'{n}\n')
        else:
            print(f"PASS: Validate {table_name}.{col_name}.")
        return num_del
