#! /usr/bin/bash
# pipelines: src/pipelines/qc/chain_clashes.py,
#            src/pipelines/dynamics/ (QC screening of geometry)

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# dssp secondary structure assignment per chain
parallel -j8 ${APP} chain_dssp ::: {0..258}

# steric / MolProbity vdW clash metrics of chains in chain_pdb
parallel -j8 ${APP} chain_clashes ::: {0..258}

# qc of original structures
parallel -j8 ${APP} pdb_qc ::: {0..258}

# qc of relaxed structures (needs 3_relaxed.sh done)
parallel -j8 ${APP} pdb_relaxed_qc ::: {0..258}
