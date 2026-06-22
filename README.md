# AlphaConformers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/AlphaConformers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/AlphaConformers.jl/dev/)
[![Build Status](https://github.com/diegozea/AlphaConformers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diegozea/AlphaConformers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/diegozea/AlphaConformers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/diegozea/AlphaConformers.jl)
[![Coverage](https://coveralls.io/repos/github/diegozea/AlphaConformers.jl/badge.svg?branch=main)](https://coveralls.io/github/diegozea/AlphaConformers.jl?branch=main)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

AlphaConformers is a Julia package for preparing structure-guided protein
prediction runs. It starts from one protein structure, searches for similar
structures, adds other known structures from related proteins, groups them by
shape, and writes prediction-ready input folders.

The goal is to help prediction tools explore alternative protein conformations,
instead of only producing models close to one preferred state.

## What It Does

AlphaConformers builds inputs for structure prediction tools in several steps:

1. It runs Foldseek on a query structure against one or more local structure
   databases.
2. It merges Foldseek hit tables and structure-based alignments.
3. It uses SIFTS and UniProt mappings to add other known structures for related
   proteins.
4. It clusters the collected templates by structural similarity.
5. It writes one folder per cluster, with an alignment and template structures.
6. It can optionally run ColabFold, AlphaFold3, or Boltz2 through helper
   functions.

The main entry point is `alphaconformers`.

## Requirements

You need:

- Julia 1.10 or newer.
- One or more local Foldseek databases.
- A folder containing PDB or mmCIF files used when adding known related
  structures.
- SIFTS mapping files, or network access so the package can download them.

If you want to run predictions from Julia, you also need the corresponding
Apptainer image and cache/model/database folders for the predictor you choose.
The input-preparation step does not require a GPU.

Useful environment variables:

- `FOLDSEEK_DB_PATH`: default Foldseek database path for low-level Foldseek
  helpers.
- `SIFTS_DB`: folder containing `pdb_chain_uniprot.csv.gz`.
- `PDB_DB`: local folder used by helper functions that read PDB files.

You can also pass paths directly to the functions instead of using environment
variables.

## Installation

From the Julia package manager:

```julia
using Pkg
Pkg.add(url = "https://github.com/diegozea/AlphaConformers.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/AlphaConformers")
```

Then load the package:

```julia
using AlphaConformers
```

## Quick Start

Prepare AlphaConformers input folders from a query structure:

```julia
using AlphaConformers

input_pdb = "1ABC_A.pdb"
pdb_folder = "datasets/pdb/mmcif_files"
out_folder = "outputs/1ABC_A"

foldseek_dbs = [
    "datasets/foldseek/fullpdb",
    "datasets/foldseek/afdb",
]

alphaconformers(
    input_pdb,
    pdb_folder,
    out_folder;
    db = foldseek_dbs,
    n_threads = 16,
    evalue_cutoff = 1e-5,
    cutoff = 1.0,
)
```

The query file name should include the PDB code and chain, for example
`1ABC_A.pdb`.

After this step, `out_folder` contains the cluster folders that can be used as
inputs for prediction.

## Running Predictions

To run ColabFold on all generated cluster folders:

```julia
using AlphaConformers

clusters_folder = "outputs/1ABC_A"
sif_path = "containers/ColabFold_AF2_1-5-5.sif"
cache_dir = "cache/colabfold"

run_alphafold(clusters_folder, sif_path, cache_dir)
```

The `cache_dir` folder is important for ColabFold runs. It is mounted inside the
Apptainer container as `/cache` and is used to store downloaded
ColabFold/AlphaFold model weights. Choose a permanent, writable location and
reuse the same folder across runs so the weights do not need to be downloaded
again. This saves time and avoids repeated large downloads. The cache can require
around 4 GB of disk space. Prediction outputs are not written to this folder; they are 
written inside each cluster's `af/` directory.

AlphaConformers also provides helpers for AlphaFold3 and Boltz2:

```julia
run_alphafold3(clusters_folder, sif_image_dir, model_parameters_dir, db_dir)
run_boltz2(clusters_folder, boltz_sif_path, boltz_cache_dir)
```

These helpers expect local container images and local model or database paths.
They do not download model weights.

## Output

The preparation step writes a folder for each template cluster:

```text
out_folder/
|-- all_sequence.a3m
|-- calpha_template.csv
|-- fullpdb_results/
|-- target_db_results/
`-- cluster_1/
    |-- sequences.a3m
    |-- templates_complete/
    `-- templates_adaptative/
```

Each `cluster_*` folder contains:

- `sequences.a3m`: the alignment used for that cluster.
- `templates_complete/`: downloaded or generated template structure files.
- `templates_adaptative/`: template files renamed for downstream prediction
  tools.

When prediction helpers are used, their raw outputs are reorganized into a
common layout:

```text
cluster_1/
`-- af/
    `-- predictions/
        `-- sequences/
            |-- models/
            |-- scores/
            `-- plots/
```

The exact prediction subfolder depends on the runner:

- `af/` for ColabFold.
- `af3/` for AlphaFold3.
- `bz/` for Boltz2.

## Common Helpers

Some lower-level helpers are useful when building custom workflows:

- `run_foldseek`: run Foldseek and save hit tables, alignments, and aligned
  structures.
- `merge_tables`: merge Foldseek hit tables.
- `merge_msas`: merge Foldseek alignments.
- `get_uniprot_mapping`: read or download SIFTS UniProt mappings.
- `create_template_clusters_hobohm`: cluster templates by structural
  similarity.
- `found_best_prediction`: filter predictions using known related structures
  and RMSD to the query.

Most users should start with `alphaconformers`.

## Development

Run the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Build the documentation:

```bash
julia --project=docs docs/make.jl
```

Source code is in `src/`, tests are in `test/`, and documentation sources are
in `docs/src/`.
