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
structural similarity, and writes prediction-ready input folders.

The goal is to help prediction tools explore alternative protein conformations,
instead of only producing models close to one preferred state.

## Installation

Download Julia from https://julialang.org/downloads/

From the Julia package manager:

```julia
using Pkg
Pkg.add(url = "https://github.com/diegozea/AlphaConformers.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/path/to/local/AlphaConformers")
```
--> path to where you have add the github link from previous step 

Then load the package:

```julia
using AlphaConformers
```
--> After that you can directly use any function mention in /src/AlphaConformers.jl

## Requirements

You need:

- Julia 1.10 or newer.
- One or more local Foldseek databases --> path save in `FOLDSEEK_DB_PATH` 
- A folder containing PDB or mmCIF files used when adding known related
  structures --> path save in `PDB_DB`
- SIFTS mapping files, or network access so the package can download them --> path save in `SIFTS_DB`

(You can also pass paths directly to the functions instead of using environment
variables)

## What It Does

AlphaConformers builds inputs for structure prediction tools in several steps:

1. It runs Foldseek on a query structure against one or more local structure
   databases --> save in `FOLDSEEK_DB_PATH` (if PDB database need to mention "pdb" in the name)
2. It merges Foldseek hit tables and structure-based alignments.
3. It uses SIFTS and UniProt mappings to add other known structures for related
   proteins --> save in `SIFTS_DB` and `PDB_DB`
4. It clusters the collected templates by structural similarity.
5. It writes one folder per cluster, with an alignment and template structures
6. It can optionally run ColabFold, AlphaFold3, or Boltz2 through helper
   functions --> path to .sif image : `COLABFOLD_SIF` and to cache directory `COLABFOLD_CACHE_DIR` where ColabFold looks for its model weights and cache files

## Quick Start

Run the default pipeline with input preparation, ColabFold prediction, and output
triage : 

```julia
using AlphaConformers

# All local variable to declare

QUERY_STRUCT = "1ABC_A.pdb"
PDB_FOLDER = "datasets/pdb/mmcif_files"
OUTPUT_DIR = "outputs/1ABC_A"

FOLDSEEK_DBS = [
    "datasets/foldseek/fullpdb",
    "datasets/foldseek/afdb",
]

COLABFOLD_SIF = "containers/colabfold-1.5.5-cuda12.2.2.sif"
COLABFOLD_CACHE_DIR = "cache/colabfold"

# Function to run the full pipeline
alphaconformers(
    COLABFOLD_SIF,
    COLABFOLD_CACHE_DIR;
    query_struct = QUERY_STRUCT,
    pdb_folder = PDB_FOLDER,
    output_dir = OUTPUT_DIR,
    databases = FOLDSEEK_DBS,
    n_threads = 16,
    evalue_cutoff = 1e-5,
    cutoff = 1.0,
)
```

The query file name should include the PDB code and chain, for example
`1ABC_A.pdb`.

## ColabFold Container Setup

AlphaConformers can run ColabFold with an Apptainer/Singularity image. A
ColabFold 1.5.5 image is available from Zenodo:

https://zenodo.org/records/20842530

Change `CONTAINER_DIR` to the folder where you want to store the image. The file
name and download URL are fixed for this Zenodo image.

```julia
using Downloads

CONTAINER_DIR = "/path/to/container/files"
COLABFOLD_SIF = joinpath(CONTAINER_DIR, "colabfold-1.5.5-cuda12.2.2.sif")

mkpath(CONTAINER_DIR)

Downloads.download(
    "https://zenodo.org/records/20842530/files/colabfold-1.5.5-cuda12.2.2.sif?download=1",
    COLABFOLD_SIF,
)
```

If your system uses Singularity instead of Apptainer, add the runtime keyword:

```julia
container_runtime = "singularity"
```

Prediction uses the default `structure_predictor = run_alphafold`. That predictor
expects two positional arguments before the semicolon: `sif_path` and `cache_dir`.
`sif_path` is the ColabFold Apptainer/Singularity `.sif` image:

```julia
COLABFOLD_SIF = "containers/colabfold-1.5.5-cuda12.2.2.sif"
```

The `cache_dir` folder is important for ColabFold runs. It is mounted inside the
Apptainer container as `/cache/colabfold`, where ColabFold looks for its model
weights and cache files. Choose a permanent location and reuse the same folder
across runs so the weights do not need to be downloaded again. This saves time
and avoids repeated large downloads. The cache can require around 4 GB of disk
space or more. Prediction outputs are not written to this folder; they are
written inside each cluster's `af/` directory.

```julia
COLABFOLD_CACHE_DIR = "cache/colabfold"
```

AlphaConformers also provides helpers for AlphaFold3 and Boltz2:

```julia
run_alphafold3(output_dir, sif_image_dir, model_parameters_dir, db_dir)
run_boltz2(output_dir, boltz_sif_path, boltz_cache_dir)
```

These helpers expect local container images and local model or database paths.
They do not download model weights. They also accept
`container_runtime = "singularity"` when needed.

## Pipeline Steps

Each step can also be run separately with `alphaconformers`. 

- Data preparation only : 

```julia
# Prepare only. Requires query_struct, pdb_folder, output_dir, and databases.
alphaconformers(;
    query_struct = QUERY_STRUCT,
    pdb_folder = PDB_FOLDER,
    output_dir = OUTPUT_DIR,
    databases = FOLDSEEK_DBS,
    predict = false,
    triage = false,
)
```
- Prediction can be run by itself:

```julia
# Predict only. Requires output_dir plus the predictor-specific arguments.
alphaconformers(
    COLABFOLD_SIF,
    COLABFOLD_CACHE_DIR;
    output_dir = OUTPUT_DIR,
    prepare = false,
    triage = false,
)
```

The input-preparation step does not require a GPU, but ColabFold
prediction is meant to run on a CUDA-capable Linux or HPC system.

- After prediction, triage can also be run separately. It requires `output_dir`,
`query_struct`, and a SIFTS mapping:

```julia
alphaconformers(;
    output_dir = OUTPUT_DIR,
    query_struct = QUERY_STRUCT,
    prepare = false,
    predict = false,
    triage = true,
    sifts_uniprot_mapping,
    foldseek_results_folder = joinpath(OUTPUT_DIR, "fullpdb_results"),
)
```
--> If `FOLDSEEK_DB_PATH` uses several databases, AlphaConformers chooses the unique database whose name contains
`pdb`, ignoring case. If it cannot choose one, it stops before running Foldseek and asks for `foldseek_results_folder`. 

## Scoring Conformers with DeepAccNet

`run_deepaccnet` runs the DeepAccNet Apptainer image over one cluster-output system
directory and writes a score table:

```julia
using AlphaConformers

base_dir = "/data/alphaconformers/1AKZ"
sif_path = "/containers/deepaccnet.sif"

csv = run_deepaccnet(base_dir, sif_path)
```

It flattens every `cluster_*/af/predictions/sequences/models/*.pdb` into a flat
folder of symlinks named `cluster_<N>_<model>.pdb`, runs the `.sif` image on a GPU
(`apptainer run --nv`), and returns the path to `deepaccnet_results.csv` (columns
`sample`, `cb-lddt`). The run needs a GPU and a working `apptainer`. Pass `dry_run = true` to print the
symlink plan and the exact apptainer command without consuming GPU time or creating
anything. Heavy intermediates (`bert_*.npy`, `*.features.npz`, `*.fa`) are removed after
a successful run unless `keep_intermediates = true`.

## Scoring Conformers with DeepAccNet

`run_deepaccnet` runs the DeepAccNet Apptainer image over one cluster-output system
directory and writes a score table:

```julia
using AlphaConformers

base_dir = "/data/alphaconformers/1AKZ"
sif_path = "/containers/deepaccnet.sif"

csv = run_deepaccnet(base_dir, sif_path)
```

It flattens every `cluster_*/af/predictions/sequences/models/*.pdb` into a flat
folder of symlinks named `cluster_<N>_<model>.pdb`, runs the `.sif` image on a GPU
(`apptainer run --nv`), and returns the path to `deepaccnet_results.csv` (columns
`sample`, `cb-lddt`). The run needs a GPU and a working `apptainer`. Pass `dry_run = true` to print the
symlink plan and the exact apptainer command without consuming GPU time or creating
anything. Heavy intermediates (`bert_*.npy`, `*.features.npz`, `*.fa`) are removed after
a successful run unless `keep_intermediates = true`.

## Output

The preparation step writes a folder for each template cluster:

```text
output_dir/
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

Full runs can take time because preparation writes one `cluster_*` folder per
structural cluster, and the default `run_alphafold` predictor runs ColabFold once
for each cluster. With the current default ColabFold settings, each processed
cluster can produce up to 10 model predictions before triage. Choose an
`output_dir` with a couple of GB free to store all those predictions. This is 
separate from `COLABFOLD_CACHE_DIR`, which stores reusable ColabFold files and
can require around 4 GB or more.

## Common Helpers

Some lower-level helpers are useful when building custom workflows:

- `run_foldseek`: run Foldseek and save hit tables, alignments, and aligned
  structures.
- `merge_tables`: merge Foldseek hit tables.
- `merge_msas`: merge Foldseek alignments.
- `get_uniprot_mapping`: read or download SIFTS UniProt mappings.
- `prepare_inputs` :  returns a small PreparedInputs object with the query path, the output folder, and the Foldseek result folder. 
- `create_template_clusters_hobohm`: cluster templates by structural
  similarity.
- `triage_outputs`: official final pipeline step for filtering and triaging
  prediction outputs.
- `found_best_prediction`: current lower-level prediction filter used by
  `triage_outputs`.

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
