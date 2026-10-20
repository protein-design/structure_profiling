#! /usr/bin/bash
# pipelines: src/pipelines/qc/chain_clashes.py,
#            src/pipelines/qc/chain_breaks.py,
#            src/pipelines/qc/chain_torsions.py,
#            src/pipelines/qc/chain_bonds.py,
#            src/pipelines/qc/chain_cbeta.py,
#            src/pipelines/qc/chain_chirality.py,
#            src/pipelines/dynamics/ (QC screening of geometry)

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# steric / MolProbity vdW clash metrics of chains in chain_pdb
parallel -j8 ${APP} chain_clashes ::: {0..258}

parallel -j8 ${APP} chain_breaks ::: {0..258}

# backbone torsion / bond-angle QC of chains in chain_pdb
parallel -j8 ${APP} chain_torsions ::: {0..258}

# covalent bond-length QC of chains in chain_pdb
parallel -j8 ${APP} chain_bonds ::: {0..258}

# MolProbity-style Cbeta-deviation QC of chains in chain_pdb
parallel -j8 ${APP} chain_cbeta ::: {0..258}

# CA-chirality QC of chains in chain_pdb
parallel -j8 ${APP} chain_chirality ::: {0..258}

# qc of original structures
# parallel -j8 ${APP} pdb_qc ::: {0..258}

# qc of relaxed structures (needs 3_relaxed.sh done)
# parallel -j8 ${APP} pdb_relaxed_qc ::: {0..258}
