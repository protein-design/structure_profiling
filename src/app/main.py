#! /usr/bin/python
import os
import sys
import subprocess

pool = {
    # dssp
    'chain_dssp': 'pdb.chain_dssp',

    # relax
    'chain_relaxed': 'dynamics.chain_relaxed',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo3_relaxed': 'dynamics.combo3_relaxed',
    'combo4_relaxed': 'dynamics.combo4_relaxed',
    'combo5_relaxed': 'dynamics.combo5_relaxed',

    # qc
    'pdb_qc': 'dynamics.pdb_qc',
    'pdb_relaxed_qc': 'dynamics.pdb_relaxed_qc',
}

def main():
    # build cmd
    name = sys.argv[1] if len(sys.argv) > 1 else None

    cmd = ['python', '-m',]
    path_prefix = 'src.pipelines.'
    if name in pool:
        cmd.append(path_prefix + pool[name])
    else:
        print("\nError: The argument <module_name> is not defined.\n")
        print("Usage: sp <module_name> [pdb_group]\n")
        print("Available pipelines:")
        for key in pool:
            print(f"  {key}")
        sys.exit(1)

    # pass arguments
    args = sys.argv[2:] if len(sys.argv) > 2 else []
    cmd += args
    print(f"Run python script: {cmd}")

    # run
    result = subprocess.run(cmd)

if __name__ == "__main__":
    main()
