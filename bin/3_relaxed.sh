#! /usr/bin/bash
# pipelines: src/pipelines/dynamics/ (relax fixed structures by OpenMM, GPU recommended)

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# one chain
parallel -j8 ${APP} chain_relaxed ::: {0..241}

# dimer
parallel -j8 ${APP} combo2_relaxed ::: {0..241}

# trimer
parallel -j8 ${APP} combo3_relaxed ::: {0..241}

# tetramer
parallel -j8 ${APP} combo4_relaxed ::: {0..241}

# pentamer
parallel -j8 ${APP} combo5_relaxed ::: {0..241}
