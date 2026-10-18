from copy import deepcopy
import mysql.connector

from .db import Db

class DeleteDb(Db):

    def __init__(self, verbose:bool=True):
        super().__init__(verbose)

    def drop_table(self, table_name):
        if self.verbose:
            print(f"Try to drop table {table_name}")
        query = f"DROP TABLE {table_name};"
        return self.execute(query)

    def empty_table(self, table_name):
        if self.verbose:
            print(f"Try to delete all records from table {table_name}...")
        query = f"DELETE FROM {table_name};"
        return self.execute(query)

    def drop_database(self, db_name=None):
        db_name = db_name if db_name else self.config.get('database')
        if not db_name:
            return False
        query = f"DROP DATABASE {db_name};"
        if self.verbose:
            print(f"Try to drop database {db_name}...")
        # connect without selecting the database being dropped
        database = self.config.pop('database', None)
        try:
            return self.execute(query)
        finally:
            if database:
                self.config['database'] = database

    def delete_rows(self, table_name, col_name:str, values:list):
        if self.verbose:
            print(f"Try to delete rows ={values} from {table_name}-{col_name}")
        placeholders = ', '.join(['%s',] * len(values))
        query = f"""
            DELETE FROM {table_name}
            WHERE {col_name} IN ({placeholders})
        """
        return self.execute(query, values)

    def delete_null_rows(self, table_name, col_name):
        query = f"""
            DELETE FROM {table_name}
            WHERE {col_name} IS NULL
        """
        return self.execute(query)