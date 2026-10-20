#! /usr/bin/bash

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# dssp secondary structure assignment per chain
parallel -j8 ${APP} chain_dssp ::: {0..258}