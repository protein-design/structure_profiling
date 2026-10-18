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
