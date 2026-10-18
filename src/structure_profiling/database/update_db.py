from copy import deepcopy
import mysql.connector

from .db import Db

class UpdateDb(Db):

    def __init__(self, verbose:bool=True):
        super().__init__(verbose)

    def execute_update(self, query, record):
        '''
        execute one query at a time
        '''
        record = deepcopy(record)
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query, record)
            connection.commit()
            return cursor.rowcount
        except mysql.connector.Error as e:
            # Roll back if the connection was established
            if 'connection' in locals() and connection.is_connected():
                connection.rollback()
            if self.verbose:
                print(f"Error: {e}")
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'cursor' in locals() and cursor is not None:
                cursor.close()
            if 'connection' in locals() and connection.is_connected():
                connection.close()
        return None

