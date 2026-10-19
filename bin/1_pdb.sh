#! /usr/bin/bash
# pipelines: src/pipelines/pdb/

APP="uv --directory /home/yuan/bio/structure_profiling run sp"

# need label_complex, meta_pdb, pdb_chain_seq, chain_pdb (built by bio_omics)
# and rawdata_dir pointing at the shared dir of downloaded PDB files

# dssp secondary structure assignment per chain
parallel -j8 ${APP} chain_dssp ::: {0..241}

# dimer
parallel -j8 ${APP} combo2 ::: {0..241}

# trimer
parallel -j8 ${APP} combo3 ::: {0..241}

# tetramer
parallel -j8 ${APP} combo4 ::: {0..241}

# pentamer
parallel -j8 ${APP} combo5 ::: {0..241}
