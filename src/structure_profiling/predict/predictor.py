'''
Execute third-party tools
Decoupled from bio_omics (bioomics.simulation.predictor):
only the colabfold runners used by the predict pipelines are kept.
'''
import os
import asyncio
from pathlib import Path
import subprocess


class Predictor:
    def __init__(self, input, outdir:Path, tools_dir:Path=None, verbose:bool=False):
        self.input = input
        self.outdir = outdir
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.tools_dir = tools_dir
        self.tool_name = None
        self.verbose = verbose
    
    @staticmethod
    async def run_command(cmd, tool_name):
        print(f"### {tool_name} asynchronous running: {cmd}")
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await process.communicate()
        returncode = process.returncode
        if returncode != 0:
            print(stderr.decode())
            return False
        print(stdout.decode())
        return True
    
    @staticmethod
    def colabfold_subprocess(params:dict, cmd:list) -> int:
        '''
        params = {
            colabfold_bin_dir: ...,
        }
        return statuscode: 0 succeed, 1-255 failed
        '''
        # conda env colabfold is needed
        env = os.environ.copy()
        env["PATH"] = (
            f"{params['colabfold_bin_dir']}:"
            + env["PATH"]
        )
        # working directory
        cwd = str(params['tools_dir'] / 'ColabFold')

        print(f"\n\nLoad {cmd[0]} command: {' '.join(cmd)}")
        try:
            res = subprocess.run(
                cmd,
                env=env,
                cwd=cwd,
                capture_output=True,
                text=True,
                check=True,
            )
            print('status: ', res.returncode)
            print('stdout: ', res.stdout)
            print('stderr: ', res.stderr)
            return res.returncode
        except subprocess.CalledProcessError as e:
            # report the real subprocess error, not only the command
            print(f"Execution colafold failed. ERROR={e}")
            print('stdout: ', e.stdout)
            print('stderr: ', e.stderr)
        except Exception as e:
            print(f"Execution colafold failed. ERROR={e}")
        return 1
    
    @staticmethod
    def colabfold_search(row:dict, params:dict):
        '''
        AlphaFold2: local MSA
        row = {
            query_faa:Path,
            outdir:Path,
        }
        '''
        cmd = [
            'colabfold_search',
            '--threads', '8',
            '--use-env', '0',
            '--use-templates', '1',
            '--db1', 'uniref30_2302_db',
            '--db2', 'pdb100_230517',
            '--db3', 'colabfold_envdb_202108_db',
        ]
        if params.get('gpu'):
            cmd += ['--gpu', '1',]
        cmd += [
            str(row['query_faa']),
            str(params['colabfold_db']),
            str(row['outdir']),
        ]        
        status = Predictor.colabfold_subprocess(params, cmd)
        if status != 0 and params.get('gpu'):
            # shared GPUs can run out of video memory (cudaMalloc/std::bad_alloc);
            # retry once without --gpu (slower CPU search)
            print('Retry colabfold_search without --gpu (CPU search)...')
            cpu_cmd = list(cmd)
            if '--gpu' in cpu_cmd:
                pos = cpu_cmd.index('--gpu')
                del cpu_cmd[pos:pos + 2]
            status = Predictor.colabfold_subprocess(params, cpu_cmd)
        return status

    @staticmethod
    def colabfold_batch_msa(row:dict, params:dict):
        '''
        AlphaFold2: access public MSA server
        row = {
            query_faa: Path,
            outdir: Path,
        }
        '''
        cmd = [
            'colabfold_batch',
            '--msa-only',
            '--templates',
            str(row['query_faa']),
            str(row['outdir']),
        ]
        return Predictor.colabfold_subprocess(params, cmd)

    @staticmethod
    def colabfold_batch(row:dict, params:dict):
        '''
        AlphaFold2: predict protein structure given sequences
        row = {
            input: Path,
            outdir: Path,
        }
        '''
        cmd = [
            'colabfold_batch',
            '--num-models', '5',
            '--model-order', '1,2,3,4,5',
            '--random-seed', '42',
        ]
        # for MSA if input is fasta
        if params.get('use_templates'):
            cmd.append('--templates')
        if params.get('save_all'):
            cmd.append('--save-all')
        if params.get('relax'):
            cmd += ['--amber', '--use-gpu-relax',]
            if 'max_iterations' in params:
                cmd += [
                    '--relax-max-iterations', str(params['max_iterations']),
                ]
        cmd += [
            '--model-type', str(params.get('model_type', 'auto')),
            '--num-recycle', str(params.get('num_recycle', 3)),
            str(row['input_file']), #fasta or a3m
            str(row['outdir']),
        ]
        return Predictor.colabfold_subprocess(params, cmd)

################################################################################

    @staticmethod
    def run_subprocess(cmd:list, cwd=None, tool_name:str=None, shell:bool=False) -> bool:
        try:
            print(f"\n\nLoad {tool_name}: {' '.join(cmd)}")
            if shell:
                cmd = ' '.join(cmd)
            res = subprocess.run(cmd, cwd=cwd, shell=shell, capture_output=True, text=True)
            print('status: ', res.returncode)
            print('stdout: ', res.stdout)
            print('stderr: ', res.stderr)
            return res
        except Exception as e:
            print(f"ERROR: failed in {tool_name} execution, error={e}")
        return None