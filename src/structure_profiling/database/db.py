'''
configurations for mysql.connector.connection
'''

import os
from pathlib import Path

import mysql.connector
from dotenv import load_dotenv

# load mysql configurations from .env
# src/.env takes precedence, the repo root .env is the fallback
load_dotenv(Path(__file__).resolve().parents[2] / '.env')
load_dotenv(Path(__file__).resolve().parents[3] / '.env')

class Db:

    def __init__(self, verbose:bool):
        self.verbose = verbose
        self.config = {
            'host': os.getenv('MYSQL_HOST'),
            'user': os.getenv('MYSQL_USER'),
            'password': os.getenv('MYSQL_PASSWORD'),
            'port': int(os.getenv('MYSQL_PORT', '3306')),
            'use_pure': True,
            'buffered': True,
            'autocommit': True,    # avoid transaction issues
            'connection_timeout': 60,
            'pool_size': 1,        # Single connection for stability
        }
        # default database for create_database/drop_database
        database = os.getenv('MYSQL_DB') or os.getenv('MYSQL_DATABASE')
        if database:
            self.config['database'] = database

    def execute(self, query, values=None) -> bool:
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query, values)
            connection.commit()
            return True
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
        return False

    def execute_fetch(self, query, values=None):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query, values)
            res = cursor.fetchall()
            connection.commit()
            return res
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'cursor' in locals() and cursor is not None:
                cursor.close()
            if 'connection' in locals() and connection.is_connected():
                connection.close()
        return -1

    def execute_many(self, query, records) -> int:
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.executemany(query, records)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'cursor' in locals() and cursor is not None:
                cursor.close()
            if 'connection' in locals() and connection.is_connected():
                connection.close()
        return 0