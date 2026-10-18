# structure_profiling

protein structure profiling for references in protein design

Standalone toolkit decoupled from
[bio_omics](https://github.com/tiezhengyuan/bio_omics): secondary structure
assignment (DSSP), OpenMM relaxation of fixed structures and QC screening of
predicted/backbone geometry.

## Layout

```
sql/                     MySQL schemas (chain_dssp, chain_fixed, chain_relaxed, combo*_fixed, combo*_relaxed, pdb_qc, pdb_relaxed_qc)
src/
├── helper.py            shared pipeline utilities (params, add_params, header, footer)
├── app/main.py          `sp` CLI: runs a pipeline as a module
├── pipelines/
│   ├── pdb/chain_dssp.py             DSSP assignment per chain
│   └── dynamics/
│       ├── chain_fixed.py            fix missing atoms/side chains (pdbfixer)
│       ├── func_fixed.py             shared fix helpers (combo2..5 reuse it)
│       ├── combo2..5_fixed.py        fix chain combos
│       ├── chain_relaxed.py          relax one chain (OpenMM)
│       ├── func_relaxed.py           shared relax helpers (combo2..5 reuse it)
│       ├── combo2..5_relaxed.py      relax chain combos
│       ├── pdb_qc.py                 QC screening of structures
│       └── pdb_relaxed_qc.py         QC screening of relaxed structures
└── structure_profiling/
    ├── common/          ProcessJson, WriteTmp (rsync tmp dir)
    ├── database/        MySQL layer (Db, Insert/Query/Delete/Update + *Complex)
    ├── structure/       Pdb parser, AnalyzeDssp
    ├── geometry/        Geometry, GeometryDistance, GeometryAngle, RamaGrid
    ├── qc/              StructureQC screening
    └── relax/           timeout, OpenmmPdbfixer, OpenmmSimulation
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

Optional: `Rama8000` grid data for `RamaGrid` goes to
`src/structure_profiling/data/Rama8000/rama8000-general-grids.npz`
(the related test is skipped when the file is absent).

## Usage

Run a pipeline (same style as bio_omics `bomics`):

```sh
sp chain_dssp 0
sp chain_fixed 0
sp combo2_fixed 0
sp chain_relaxed 0
sp pdb_qc 0
sp pdb_relaxed_qc 0
```

`0` is the `pdb_group` to process; when omitted it is auto-assigned from the
`meta_execute` table.

Create the tables once per database:

```sh
mysql COMPLEX2 < sql/pdb.sql
mysql COMPLEX2 < sql/simulate.sql
```

## Tests

```sh
uv run pytest tests/ -v
```
