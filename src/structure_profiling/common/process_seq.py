'''
process sequences exported as the input of structure prediction.
Decoupled from bio_omics (bioomics.protein.process_seq): only the
colabfold exporters used by the combo2..5 pipelines are kept.
'''


class ProcessSeq:

    @staticmethod
    def to_colabfold_csv(data:dict, outfile:str):
        '''
        Used as the input of colabfold for structure prediction
        '''
        try:
            with open(outfile, 'w') as f:
                f.write("id,sequence\n")
                for seq_id, seq_pool in data.items():
                    seq_str = ':'.join(seq_pool)
                    text = f"{seq_id},{seq_str}\n"
                    f.write(text)
            return 'succeed'
        except Exception as e:
            return 'fail'

    @staticmethod
    def to_colabfold_fasta(data:dict, outfile:str):
        '''
        Used as the input of colabfold for structure prediction
        '''
        try:
            with open(outfile, 'w') as f:
                for seq_id, seq_pool in data.items():
                    header = f">{seq_id}\n"
                    f.write(header)
                    seq_str = ':\n'.join(seq_pool)
                    f.write(seq_str)
            return 'succeed'
        except Exception as e:
            return 'fail'