# structure_profiling

protein structure profiling for references in protein design

Standalone toolkit decoupled from
[bio_omics](https://github.com/tiezhengyuan/bio_omics): split of protein
complexes into chain combos, secondary structure assignment (DSSP), OpenMM
relaxation of fixed structures, AlphaFold (colabfold) structure prediction of
chains/combos and QC screening of predicted/backbone geometry.

## Layout

```
sql/                     MySQL schemas (chain_dssp, chain_fixed, chain_relaxed, combo2..5, combo*_fixed, combo*_relaxed, chain_clashes, af_params, af_msa, af_pred, af_summary, af_plddt, af_rmsd, af_tmalign)
src/
├── helper.py            shared pipeline utilities (params, add_params, header, footer)
├── app/main.py          `sp` CLI: runs a pipeline as a module
├── pipelines/
│   ├── pdb/
│   │   ├── chain_dssp.py             DSSP assignment per chain
│   │   ├── func_combo.py             shared combo helpers (combo2..5 reuse it)
│   │   └── combo2..5.py              split complexes into chain combos (csv/fasta/pdb export)
│   └── dynamics/
│       ├── chain_fixed.py            fix missing atoms/side chains (pdbfixer)
│       ├── func_fixed.py             shared fix helpers (combo2..5 reuse it)
│       ├── combo2..5_fixed.py        fix chain combos
│       ├── chain_relaxed.py          relax one chain (OpenMM)
│       ├── func_relaxed.py           shared relax helpers (combo2..5 reuse it)
│       └── combo2..5_relaxed.py      relax chain combos
│   ├── qc/
│   │   ├── chain_clashes.py           steric / MolProbity vdW clash metrics per chain
│   │   ├── chain_breaks.py            chain continuity (missing atoms, breaks, strain)
│   │   ├── chain_fixed_breaks.py      chain continuity of fixed chains (input table chain_fixed)
│   │   ├── chain_torsions.py          torsions, bond-angle Z, Ramachandran, planarity
│   │   ├── chain_bonds.py             covalent bond-length Z metrics
│   │   ├── chain_cbeta.py             Cbeta-deviation metrics
│   │   └── chain_chirality.py         CA-chirality (handedness) metrics
│   ├── predict/
│   │   ├── func_predict.py             shared prediction helpers (af_params/af_msa/af_pred reuse it)
│   │   ├── af_params.py                build prediction inputs from chain/combo sequences
│   │   ├── af_msa.py                   build MSA (*.a3m) via colabfold
│   │   ├── af_pred.py                  predict structures with colabfold batch
│   │   ├── af_summary.py               per-model summary (ranking, pLDDT, PAE)
│   │   ├── af_plddt.py                 per-chain pLDDT of predicted models
│   │   ├── af_rmsd.py                  RMSD between predictions and native structures
│   │   ├── af_tmalign.py               TM-align score between predictions and native structures
│   │   └── func_summary.py             shared summary helpers (af_rmsd/af_tmalign reuse it)
└── structure_profiling/
    ├── common/          ProcessJson, ProcessSeq (colabfold csv/fasta export), WriteTmp (rsync tmp dir)
    ├── database/        MySQL layer (Db, Insert/Query/Delete/Update + *Complex)
    ├── structure/       Pdb parser (pdb.py, parse_structure, atom_select;
    │                    to_chain_pdb), AnalyzeDssp, timeout, OpenmmPdbfixer,
    │                    OpenmmSimulation
    └── geometry/        Breaks, Torsions, Clashes + protein-geometry QC:
                         Geometry, GeometryDistance, GeometryAngle, RamaGrid
tests/
```

## Install

```sh
uv sync          # creates .venv, installs deps (pdbfixer from git)
cp .env.example .env   # then fill in MySQL credentials and output_dir
```

For GPU hosts install the CUDA build of OpenMM instead: `uv add "openmm[cuda12]"`.

DSSP assignment additionally requires the `mkdssp` executable (e.g.
`apt install dssp` or `conda install -c salilab dssp`).

The `af_msa` / `af_pred` predict pipelines additionally require a colabfold
installation (bin dir and databases configured in `.env`, see `.env.example`).

Optional: `Rama8000` grid data for `RamaGrid` goes to
`src/structure_profiling/data/Rama8000/rama8000-general-grids.npz`
(the related test is skipped when the file is absent).

## Usage

Run a pipeline (same style as bio_omics `bomics`):

```sh
sp test_env
sp chain_dssp 0
sp combo2 0
sp chain_fixed 0
sp combo2_fixed 0
sp chain_relaxed 0
sp chain_clashes 0
sp chain_breaks 0
sp chain_fixed_breaks 0
sp af_params 0
sp af_msa 0
sp af_pred 0
sp af_summary 0
sp af_plddt 0
sp af_rmsd 0
sp af_tmalign 0
```

`sp test_env` takes no group: it reports python and dependency versions, the
`.env` in use (secrets masked, output dirs checked), OpenMM platforms and
CUDA/GPU availability (PyTorch when installed), `mkdssp` plus its runtime CCD
data (`components.cif`, which mkdssp 4.x needs to resolve residue compounds),
and MySQL connectivity. It exits non-zero when a required dependency or path
is broken.

`0` is the `pdb_group` to process; when omitted it is auto-assigned from the
`meta_execute` table.

The `combo2..5` split pipelines need `label_complex`, `meta_pdb`, `pdb_chain_seq`
and `chain_pdb` filled by bio_omics, plus the `rawdata_dir` env var pointing at
the shared directory of downloaded PDB files (`meta_pdb.relative_pdb`). Their
output feeds the `combo2..5_fixed` / `combo2..5_relaxed` pipelines.

The `af_*` predict pipelines need `label_chain` / `label_complex` and the
`pdb_chain_seq` fasta files (both filled by bio_omics, which also creates the
shared `pdb` table referenced by `af_params`). MSA/prediction require a
colabfold installation (`colabfold_bin_dir`, `colabfold_db` in `.env`),
`af_tmalign` needs TMalign (`mmalign_dir`); predictions are written below
`predict_dir` (default `output_dir/predict`).

Create the tables once per database:

```sh
mysql COMPLEX2 < sql/pdb.sql
mysql COMPLEX2 < sql/simulate.sql
mysql COMPLEX2 < sql/qc.sql
mysql COMPLEX2 < sql/predict.sql
```

## Tests

```sh
uv run pytest tests/ -v
```
