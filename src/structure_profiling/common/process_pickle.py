'''
binary file
Decoupled from bio_omics (bioomics.common.process_pickle).
'''
from pathlib import Path
import pickle

class ProcessPickle:

    def __init__(self, pickle_file, verbose:bool=False):
        # NB: coerce to Path so that str inputs work with .is_file()
        self.pickle_file = Path(pickle_file)
        self.verbose = verbose
    
    def save_data(self, data:dict):
        if self.pickle_file.is_file():
            old_data = self.load_pickle()
            old_data.update(data)
            return self.save_pickle(old_data)
        return self.save_pickle(data)

    def save_pickle(self, data) -> str:
        try:
            with open(self.pickle_file, 'wb') as f:
                pickle.dump(data, f)
            if self.verbose:
                print(f"Save data to {self.pickle_file}")
            return 'succeed'
        except Exception as e:
            if self.verbose:
                print(f"Error: Failed to save data to pickle {self.pickle_file}, error={e}")
            pass
        return 'fail'
        return 'fail'

    def load_pickle(self, default=None):
        try:
            with open(self.pickle_file, 'rb') as f:
                if self.verbose:
                    print(f"Try to read data from {self.pickle_file}.")
                return pickle.load(f)
        except Exception as e:
            if self.verbose:
                print(f"Warning: No data loading from pickle {self.pickle_file}, error={e}")
            pass
        return default
    
    def is_good(self)-> bool:
        '''
        return False if empty, not existing, or wrong pickle file
        '''
        data = self.load_pickle()
        if data:
            return True
        return False