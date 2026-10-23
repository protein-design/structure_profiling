/*
chain_ligand: ligand / binding-site neighborhood per chain
(pipelines/dynamics/chain_ligand.py).

For every chain of chain_pdb the HETATM ligands (waters skipped) are
taken as sites - or an explicit activity site given by site residues -
and every standard residue of the chain with a heavy atom within
search_cutoff (default 5.0 A, structure/neighbors.py) is a neighbor.
Neighbors may be scattered over the sequence (pocket): n_fragments
counts the sequence-continuous segments among them. The full contact
table goes to neighbors.csv / neighbors.pkl below output_pdb_dir
(neighbor_csv / neighbor_pkl hold the relative paths).

Schema references the shared COMPLEX2 tables (pdb, chain_pdb).
*/

DROP TABLE IF EXISTS chain_ligand;
CREATE TABLE chain_ligand(
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) PRIMARY KEY,
    n_sites             INT,
    ligand_names        VARCHAR(200),
    n_neighbor_residues INT,
    n_fragments         INT,
    min_dist            FLOAT,
    mean_dist           FLOAT,
    neighbor_residues   TEXT,
    neighbor_csv        VARCHAR(120),
    neighbor_pkl        VARCHAR(120),
    CONSTRAINT fk_chainligand_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);