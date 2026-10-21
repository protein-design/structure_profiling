'''
read FASTA files
Decoupled from bio_omics (bioomics.protein.process_fasta):
only the reader used by the predict pipelines is kept
(writers stay in ProcessSeq).
'''
from Bio import SeqIO
from pathlib import Path

class ProcessFasta:
    def __init__(self, fa_file:Path, verbose:bool=False):
        self.fa_file = Path(fa_file)
        self.verbose = verbose

    def load_fasta(self) -> list:
        try:
            with open(self.fa_file, 'r') as f:
                parser = SeqIO.parse(f, 'fasta')
                records = [rec for rec in parser]
                return records
        except Exception as e:
            if self.verbose:
                print(f"Failure: Can't read FASTA {self.fa_file}, error={e}")
        return []

    def get_first_seq(self) -> str:
        for rec in self.load_fasta():
            return str(rec.seq)
        return ''