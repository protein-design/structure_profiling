#! /usr/bin/python
import os
import sys
import subprocess

pool = {
    # diagnostics
    'test_env': 'test_env',

    # pdb
    'chain_dssp': 'pdb.chain_dssp',

    # pdb: split complex into chain combos
    'combo2': 'pdb.combo2',
    'combo3': 'pdb.combo3',
    'combo4': 'pdb.combo4',
    'combo5': 'pdb.combo5',

    # fix
    'chain_fixed': 'dynamics.chain_fixed',
    'combo2_fixed': 'dynamics.combo2_fixed',
    'combo3_fixed': 'dynamics.combo3_fixed',
    'combo4_fixed': 'dynamics.combo4_fixed',
    'combo5_fixed': 'dynamics.combo5_fixed',

    # relax
    'chain_relaxed': 'dynamics.chain_relaxed',
    'combo2_relaxed': 'dynamics.combo2_relaxed',
    'combo3_relaxed': 'dynamics.combo3_relaxed',
    'combo4_relaxed': 'dynamics.combo4_relaxed',
    'combo5_relaxed': 'dynamics.combo5_relaxed',

    # qc
    'chain_clashes': 'qc.chain_clashes',
    'chain_breaks': 'qc.chain_breaks',
    'chain_torsions': 'qc.chain_torsions',
    'chain_bonds': 'qc.chain_bonds',
    'chain_cbeta': 'qc.chain_cbeta',
    'chain_chirality': 'qc.chain_chirality',
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
