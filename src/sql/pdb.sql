/*
chain_dssp: DSSP secondary structure assignment per chain
Moved from bio_omics src/pipelines/create/pdb.sql.

Schema references the shared COMPLEX2.pdb table.
*/

DROP TABLE IF EXISTS chain_dssp;
CREATE TABLE chain_dssp(
    pdb_id    VARCHAR(10),
    chain_id  VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    model_no  INT,
    chain_no  VARCHAR(5),
    dssp_seq  TEXT,
    host      VARCHAR(10),
    UNIQUE(chain_id),
    CONSTRAINT fk_chaindssp_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
