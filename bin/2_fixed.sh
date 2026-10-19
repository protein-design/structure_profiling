#! /usr/bin/bash
# pipelines: src/pipelines/dynamics/ (fix missing atoms/side chains by pdbfixer)

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# one chain
parallel -j8 ${APP} chain_fixed ::: {0..241}

# dimer
parallel -j8 ${APP} combo2_fixed ::: {0..241}

# trimer
parallel -j8 ${APP} combo3_fixed ::: {0..241}

# tetramer
parallel -j8 ${APP} combo4_fixed ::: {0..241}

# pentamer
parallel -j8 ${APP} combo5_fixed ::: {0..241}
