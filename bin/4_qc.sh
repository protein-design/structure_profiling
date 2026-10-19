#! /usr/bin/bash
# pipelines: src/pipelines/dynamics/ (QC screening of geometry)

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# qc of original structures
parallel -j8 ${APP} pdb_qc ::: {0..241}

# qc of relaxed structures (needs 3_relaxed.sh done)
parallel -j8 ${APP} pdb_relaxed_qc ::: {0..241}
