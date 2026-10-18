'''
basic operation for altering tables
'''
import itertools
import mysql.connector

from .db import Db

class InsertDb(Db):

    def __init__(self, verbose:bool=False, chunk_size=None):
        Db.__init__(self, verbose)
        self.chunk_size = int(chunk_size) if chunk_size else 100
   
    def insert(self, data:list, table_name:str, cols:list):
        '''
        data is one record
        '''
        col_str = ', '.join(cols)
        fmt_str = ', '.join(['%s'] * len(data))
        query = f"""
            INSERT INTO {table_name} ({col_str}) VALUES ({fmt_str})
        """
        return self.execute(query, tuple(data))

    def insert_records(self, data:list, table_name:str, cols:list) -> tuple:
        '''
        data are records in iterable type
        '''
        cols_str = ', '.join(cols)
        fmt = ['%s',] * len(cols)
        fmt_str= ', '.join(fmt)
        query = f"""
            INSERT INTO {table_name} ({cols_str})
            VALUES ({fmt_str})
        """
        succeed, fail = 0, []
        for chunk_records in itertools.batched(data, self.chunk_size):
            counts = self.execute_many(query, chunk_records)
            if counts:
                succeed += counts
            else:
                for record in chunk_records:
                    n = self.execute_many(query, [record,])
                    if n:
                        succeed += n
                    else:
                        fail.append(record)
        return succeed, fail
    
    def func_insert(self, data:list, insert_query, func, **kwargs) -> tuple:
        '''
        data is parameters, func would build records from self.data
        return unscceessful data if some insert failed
        '''
        # get records: self.data is pdb paths
        n, succeed = 0, 0
        pool, recycle_pool = [], []
        for arg in data:
            n += 1
            if n % 1_000 == 0:
                print(n, end='..')
            # collect record: arg is pdb_file or dat_record
            records = func(arg, **kwargs)
            pool.extend(records)
            # insert records into table
            if len(pool) >= self.chunk_size:
                m = self.execute_many(insert_query, pool)
                if m:
                    succeed += len(pool)
                else:
                    recycle_pool += pool
                pool = []
        if pool:
            m = self.execute_many(insert_query, pool)
            if m:
                succeed += len(pool)
            else:
                recycle_pool += pool
        print(f"\nTry again: insert one by one...")
        failed = []
        for rec in recycle_pool:
            m = self.execute_many(insert_query, [rec,])
            if m:
                succeed += 1
            else:
                failed.append(rec)
        print(f"failed:{len(failed)}, succeed={succeed}.\n")
        return succeed, failed
