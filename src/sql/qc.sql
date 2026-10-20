/*
Chain-clash, chain-continuity and chain-torsion QC tables
(pipelines/qc/chain_clashes.py, pipelines/qc/chain_breaks.py,
pipelines/qc/chain_torsions.py).

Cutoffs live in structure/clashes.py (Clashes), qc/breaks.py (Breaks)
and qc/torsions.py (Torsions).
*/

DROP TABLE IF EXISTS chain_clashes;
CREATE TABLE chain_clashes(
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) PRIMARY KEY,
    n_heavy_atoms       INT,
    num_clashes         INT,
    num_severe_clashes  INT,
    clash_score         FLOAT,
    n_vdw_clashes       INT,
    vdw_clash_score     FLOAT,
    max_overlap         FLOAT,
    CONSTRAINT fk_chainclashes_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

DROP TABLE IF EXISTS chain_breaks;
CREATE TABLE chain_breaks(
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) PRIMARY KEY,
    n_residues          INT,
    num_missing         INT,
    num_breaks          INT,
    num_fragments       INT,
    num_strained        INT,
    ca_max_dist         FLOAT,
    pep_max_dist        FLOAT,
    CONSTRAINT fk_chainbreaks_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

/* backbone torsion / bond-angle QC of chains (pipelines/qc/chain_torsions.py).
   Cutoffs live in qc/torsions.py (Torsions):
   Z_OUTLIER_CUT=5, CIS_CUT=90 deg, OMEGA_DEV_CUT=30 deg,
   MAINCHAIN_OOP_CUT=10 deg; GATES: angle_rmsz=1.5, z_outlier_frac=0.01,
   torsion_null_frac=0.05, peptide_outlier_frac=0.05,
   mainchain_outlier_frac=0.05, rama_outlier_frac=0.05 */
DROP TABLE IF EXISTS chain_torsions;
CREATE TABLE chain_torsions(
    pdb_id                  VARCHAR(10),
    chain_id                VARCHAR(20) PRIMARY KEY,
    n_residues              INT,
    n_angles                INT,
    angle_rmsz              FLOAT,
    n_z_outliers            INT,
    z_outlier_frac          FLOAT,
    n_torsions              INT,
    torsion_null_frac       FLOAT,
    frac_trans              FLOAT,
    frac_cis                FLOAT,
    n_cis_nonpro            INT,
    rama_analysed           INT,
    rama_favoured_frac      FLOAT,
    rama_allowed_frac       FLOAT,
    rama_outlier_frac       FLOAT,
    n_peptides              INT,
    peptide_outlier_frac    FLOAT,
    mainchain_outlier_frac  FLOAT,
    CONSTRAINT fk_chaintorsions_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);



DROP TABLE IF EXISTS pdb_qc;
CREATE TABLE pdb_qc(
    pdb_id                  VARCHAR(10),
    chain_id                VARCHAR(20) PRIMARY KEY,
    angle_rmsz              FLOAT,
    bond_rmsz               FLOAT,
    ca_max_dist             FLOAT,
    cb_outlier_frac         FLOAT,
    clash_score             FLOAT,
    mainchain_outlier_frac  FLOAT,
    n_chirality_outliers    INT,
    n_z_outliers            INT,
    num_breaks              INT,
    num_clashes             INT,
    num_fragments           INT,
    num_missing             INT,
    num_strained            INT,
    num_servere_clashes     INT,
    pep_max_dist            FLOAT,
    peptide_outlier_frac    FLOAT,
    rama_outlier_frac       FLOAT,
    vdw_clash_score         FLOAT,
    z_outlier_frac          FLOAT,
    CONSTRAINT fk_pdbqc_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);
CREATE TABLE pdb_relaxed_qc LIKE pdb_qc;

