import gc
import pandas as pd
import time

from .db import Db

class QueryDb(Db):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)

    def list_data(self, query:str, is_df:bool=False):
        # try not more than 3 times
        for _ in range(3):
            res = self.execute_fetch(query)
            if res == -1:
                gc.collect()
                time.sleep(5)
            else:
                if is_df:
                    return pd.DataFrame(res)
                return res
        return []

    def parse_data(self, query:str, values:tuple) -> list:
        placeholders = ', '.join(['%s'] * len(values))
        query += f"({placeholders})"
        # try not more than 3 times
        for _ in range(3):
            res = self.execute_fetch(query, tuple(values))
            if res == -1:
                gc.collect()
                time.sleep(5)
            else:
                return res
        return []
    
    def table_exists(self, table_name:str) -> bool:
        res = self.execute_fetch("SHOW TABLES")
        if not isinstance(res, list) or len(res) == 0:
            return False
        tables = pd.DataFrame(res).iloc[:,0].to_list()
        return table_name in tables

    def num_rows(self, table_name:str) -> int:
        query = f"""
            SELECT COUNT(*) AS count FROM {table_name}
        """
        rows = self.execute_fetch(query)
        if rows and rows != -1:
            return rows[0]['count']
        return -1

    def unique_values(self, table_name, col_name:str) -> list:
        query = f"""
            SELECT DISTINCT {col_name}
            FROM {table_name}
            WHERE {col_name} IS NOT NULL;
        """
        rows = self.list_data(query)
        if self.verbose:
            print(f"Detect {len(rows)} records in {table_name}.{col_name}.")
        return [i[col_name] for i in rows]

    def next_id(self, table_name, col_name):
        query = f"""
            SELECT MAX({col_name}) AS id FROM {table_name};
        """
        res = self.list_data(query)
        curr = res[0]['id']
        return int(curr) + 1 if curr else 1

    def max_value(self, table_name:str, col:str):
        query = f"""
            SELECT MAX({col}) AS val
            FROM {table_name}
        ;"""
        res = self.list_data(query)
        return res[0]['val']

    def table_names(self):
        query = "SHOW TABLES;"
        res = self.list_data(query)
        if not isinstance(res, list) or len(res) == 0:
            return ()
        res = pd.DataFrame(res).iloc[:,0]
        return tuple(res)
