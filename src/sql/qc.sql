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

/* covalent bond-length QC of chains (pipelines/qc/chain_bonds.py).
   Cutoffs live in geometry/geometry_distance.py (GeometryDistance):
   Z_OUTLIER_CUT=5 on |Z| vs the Engh & Huber IDEAL_BONDS;
   former BACKBONE_GATES: bond_rmsz=2.0 (wwPDB expects [0, 1]),
   z_outlier_frac=0.01 (bond + angle rows combined) */
DROP TABLE IF EXISTS chain_bonds;
CREATE TABLE chain_bonds(
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) PRIMARY KEY,
    n_bonds             INT,
    bond_rmsz           FLOAT,
    n_bond_outliers     INT,
    bond_z_outlier_frac FLOAT,
    bond_max_absz       FLOAT,
    CONSTRAINT fk_chainbonds_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

/* MolProbity-style Cbeta-deviation QC of chains
   (pipelines/qc/chain_cbeta.py). Cutoffs live in
   geometry/geometry_distance.py (GeometryDistance):
   CB_DEVIATION_CUT=0.25 A; former STRUCTURE_GATES cb_outlier_frac=0.05 */
DROP TABLE IF EXISTS chain_cbeta;
CREATE TABLE chain_cbeta(
    pdb_id              VARCHAR(10),
    chain_id            VARCHAR(20) PRIMARY KEY,
    n_residues          INT,
    n_cb_analysed       INT,
    cb_max_dev          FLOAT,
    cb_outlier_frac     FLOAT,
    CONSTRAINT fk_chaincbeta_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);

/* CA-chirality QC of chains (pipelines/qc/chain_chirality.py).
   Cutoffs live in geometry/geometry_angle.py (GeometryAngle):
   CHIRAL_VOL_L=2.6 A^3, CHIRAL_VOL_MIN=0.7 A^3 (0 < volume < min is
   'planar', < 0 is D, no CB is 'NA'); former STRUCTURE_GATES
   n_chirality_outliers=0 */
DROP TABLE IF EXISTS chain_chirality;
CREATE TABLE chain_chirality(
    pdb_id                  VARCHAR(10),
    chain_id                VARCHAR(20) PRIMARY KEY,
    n_residues              INT,
    n_chirality_analysed    INT,
    n_chirality_outliers    INT,
    chirality_outlier_frac  FLOAT,
    CONSTRAINT fk_chainchirality_pdb FOREIGN KEY (pdb_id)
        REFERENCES pdb(pdb_id) ON DELETE CASCADE
);



