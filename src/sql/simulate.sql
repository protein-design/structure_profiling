/*
Fixed-, relaxed-structure and QC tables.
Moved from bio_omics src/pipelines/create/simulate.sql.

meta_pdb -> chain_pdb -> chain_fixed -> chain_relaxed
         -> combo2 -> combo2_fixed -> combo2_relaxed
         -> combo3 -> combo3_fixed -> combo3_relaxed
         -> combo4 -> combo4_fixed -> combo4_relaxed
         -> combo5 -> combo5_fixed -> combo5_relaxed
         -> chain_clashes -> chain_breaks -> chain_torsions
         -> chain_bonds -> chain_cbeta -> chain_chirality

Schema references the shared COMPLEX2 tables (pdb, combo2..combo5).
*/

DROP TABLE IF EXISTS chain_fixed;
CREATE TABLE chain_fixed (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_relaxed;
CREATE TABLE chain_relaxed (
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) COLLATE utf8mb4_bin NOT NULL,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (chain_id),
    CONSTRAINT fk_chain_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo2_fixed;
CREATE TABLE combo2_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo2_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo2_relaxed;
CREATE TABLE combo2_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo2_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo3_fixed;
CREATE TABLE combo3_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo3_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo3_relaxed;
CREATE TABLE combo3_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo3_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo4_fixed;
CREATE TABLE combo4_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo4_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo4_relaxed;
CREATE TABLE combo4_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo4_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo5_fixed;
CREATE TABLE combo5_fixed(
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    fixed_pdb           VARCHAR(80) NOT NULL,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo5_fixed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS combo5_relaxed;
CREATE TABLE combo5_relaxed (
    pdb_id              VARCHAR(10),
    combo_id            VARCHAR(40) COLLATE utf8mb4_bin NOT NULL,
    relaxed_pdb         VARCHAR(80),
    original_energy     VARCHAR(60),
    relaxed_energy      VARCHAR(60),
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE (combo_id),
    CONSTRAINT fk_combo5_relaxed_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);


