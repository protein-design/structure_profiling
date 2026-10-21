'''
requirements:
    - tables af_pred
example: 
    - parallel -j8 sp af_summary ::: {0..241}
    - sp af_summary 0
functions:
    - retrieve summary of alphafold prediction to table af_summary
'''

from src.helper import *
from structure_profiling import ProcessJson
from .func_predict import *
from .func_summary import *
            
def parse_data(row):
    res = {}
    for path in Path(row['pred_outdir']).iterdir():
        if 'rank_' in path.name:
            ranking = int(re.findall(r'rank_(\d+)_', path.name)[0])
            if ranking not in res:
                res[ranking] = {'ranking': ranking,}
            if '_scores_' in path.name and path.suffix == '.json':
                data = ProcessJson.load_json(path)
                res[ranking]['avg_plddt'] = np.mean(data['plddt'])
                res[ranking]['avg_ptm'] = np.mean(data['ptm'])
                res[ranking]['max_pae'] = data['max_pae']
                res[ranking]['json_file_name'] = path.name
            elif '_relaxed_' in path.name and path.suffix == '.pdb':
                res[ranking]['relaxed_pdb_name'] = path.name
            elif '_unrelaxed_' in path.name and path.suffix == '.pdb':
                res[ranking]['unrelaxed_pdb_name'] = path.name
    return res

def build(rows, pdb_names, meta):
    for row in rows:
        ranked = parse_data(row)
        for ranking, log in ranked.items():
            relaxed_pdb_name = log.get('relaxed_pdb_name')
            if relaxed_pdb_name not in pdb_names:
                rec = (
                    row['input_id'],
                    ranking,
                    log.get('avg_plddt'),
                    log.get('avg_ptm'),
                    log.get('max_pae'),
                    log.get('json_file_name'),
                    relaxed_pdb_name, 
                    log.get('unrelaxed_pdb_name'),
                )
                print(rec)
                yield rec
                meta['records'] += len(ranked)

if __name__ == "__main__":
    params.update({
        'params_table_name': 'af_params',
        'pred_table_name': 'af_pred',
        'table_name': 'af_summary',
    })
    
    # detect predictions from af_pred
    rows = retrieve_pred(params, meta)
    
    # detect files in af_summary
    relaxed_pdb_names = pred_summary(params)

    # parse results
    record_iter = build(rows, relaxed_pdb_names, meta)

    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
    
