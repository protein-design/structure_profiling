'''
update records of the COMPLEX2 database
Subset of bio_omics bioomics.database.update_complex: only the updates
used by the structure profiling pipelines are kept.
'''
import os
from datetime import datetime
from .update_db import UpdateDb

class UpdateComplex(UpdateDb):

    def __init__(self, verbose:bool=False):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')

    def execute_end_time(self, args):
        t1 = datetime.now()
        query = f"""
            update meta_execute set end_at = '{t1}'
            where table_name = '{args['table_name']}'
                and pdb_group = {args['pdb_group']}
                and host = '{args['host']}'
        ;"""
        self.execute_update(query, None)
