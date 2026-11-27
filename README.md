# AlphaConformers

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.08.18.745512-b31b1b.svg)](https://www.biorxiv.org/content/10.64898/2026.08.18.745512v1)
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

Zenodo link for reproducibility: https://zenodo.org/records/21133448

## Summary

- [Installation](#installation)
- [Requirements](#requirements)
- [What It Does](#what-it-does)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Scoring Conformers with DeepAccNet](#scoring-conformers-with-deepaccnet)
- [Conformer Clustering](#conformer-clustering)
- [Container Setup](#container-setup)
- [Common Helpers](#common-helpers)
- [Development](#development)
- [Citation](#citation)

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
- A folder containing PDB or mmCIF files used as input for AF2 predictions --> path save in `PDB_DB`.
- SIFTS mapping files, or network access so the package can download them automatically.

## What It Does

AlphaConformers builds inputs for structure prediction tools in several steps:

1. It runs Foldseek on a query structure against one or more local structure
   databases --> save in `FOLDSEEK_DB_PATH` (if PDB database need to mention "pdb" in the name)
2. It merges Foldseek hit tables and structure-based alignments.
3. It uses SIFTS and UniProt mappings to add other known structures for related
   proteins --> save in `PDB_DB`
4. It clusters the collected templates by structural similarity.
5. It writes one folder per cluster, with an alignment and template structures
6. It can optionally run ColabFold, AlphaFold3, or Boltz2 through helper
   functions --> path to .sif image : `COLABFOLD_SIF` and to cache directory `COLABFOLD_CACHE_DIR` where ColabFold looks for its model weights and cache files

More details on each step are in the [preprint](https://www.biorxiv.org/content/10.64898/2026.08.18.745512v1).

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

Running `alphaconformers` writes results into `output_dir/`:

```text
output_dir/
|-- fullpdb_results/
|-- all_sequence.a3m
`-- cluster_1/
    |-- sequences.a3m
    |-- templates_complete/
    |-- templates_adaptative/
    `-- af/predictions/sequences/models/   # predicted structures (.pdb)
```

Each `cluster_*` folder contains:

- `sequences.a3m`: the alignment used for that cluster.
- `templates_complete/`: downloaded or generated template structure files.
- `templates_adaptative/`: template files renamed for downstream prediction tools.

Full runs can take hours because preparation creates one `cluster_*` folder per
structural cluster and each cluster can produce up to 10 predicted models. Choose
an `output_dir` with a couple of GB of free space to store all the predictions.

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

AlphaConformers also provides helpers for AlphaFold3, Boltz2.

```julia
using AlphaConformers

run_alphafold3(clusters_folder, sif_image_dir, model_parameters_dir, db_dir)
run_boltz2(clusters_folder, boltz_sif_path, boltz_cache_dir)
```

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

`run_deepaccnet` scores the predicted conformers in one AlphaConformers output
directory. It runs the DeepAccNet Apptainer/Singularity image and writes a CSV
score table (see [Container Setup](#container-setup) below to download
`deepaccnet.sif` from Zenodo).

```julia
using AlphaConformers

output_dir = "/data/alphaconformers/pdb/1AKZ"
sif_path = DEEPACCNET_SIF

csv = run_deepaccnet(output_dir, sif_path)
```

It auto-detects the inner prediction-method folder, runs the `.sif` image on a GPU
(`apptainer run --nv`), and returns the path to
`output_dir/deepaccnet/deepaccnet_results.csv`. The run needs a GPU and a working
`apptainer` (or `singularity`, via `container_runtime`). Pass `n_threads` to
control the number of CPUs used for featurization.

## Conformer Clustering

`cluster_conformers` groups an AlphaConformers prediction ensemble for one
protein system by shape and writes a per-system output tree plus the essential
figures. It runs in three stages:

1. **KMeans clustering**: groups predicted conformers by structural similarity.
2. **DeepAccNet score filter** (optional): filters clusters using a score table.
3. **Agglomerative sub-clustering**: sub-clusters each occupied (or surviving) KMeans cluster.

Reference-free clustering:

```julia
using AlphaConformers

cluster_conformers("/data/alphaconformers/1AKZ_A")
```

With one or two reference structures:

```julia
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb")
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb", "refs/1SSP_E.pdb")
```

With common options (shown with their defaults):

```julia
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb", "refs/1SSP_E.pdb";
    method     = "af",             # prediction-method subfolder to read (alphaconformers folder_af2_result)
    ref_labels = ("apo", "holo"))  # label each reference (default: file basename)
```

Clustering notes:

- **Run folder**: the single positional argument is the `output_dir` written by
  `alphaconformers`. Predictions are read from
  `cluster_*/<method>/predictions/sequences/models/`, and clustering results are
  written into `conformer_clustering/` inside that same folder.
- **References**: each reference is a path to a `.pdb`/`.cif` file (optionally
  `.gz`) with no chain argument. The reference count (0/1/2) is auto-detected,
  and labels come from `ref_labels` (default: file basename without extension) to
  name `RMSD_<label>` and the cluster-assignment row.
- **Defaults**: clustering parameters (`K = 30`, agglomerative cut `1.0 Å`,
  average linkage, seed `42`) are fixed defaults.

To run the optional DeepAccNet score filter (Stage 2), pass a score table - a
`DataFrame` or a path to a CSV with a `Name` column and the score column:

```julia
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb", "refs/1SSP_E.pdb"; score_table = "deepaccnet_results.csv")
```

### Query path

The query-path figure (`query_path.png`) ranks the clusters by their RMSD to a
**query** structure, starting from the cluster the query falls into. By default the
query is the first reference; with no reference and no `query`, the figure is skipped.

```julia
# query_path.png is produced automatically - the first reference is the query.
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb")                     # 1 reference
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb", "refs/1SSP_E.pdb")  # 2 references → first is the query
```

To override the default, pass `query` explicitly - either a conformer name (as it
appears in the clustering table) or a path to a structure file (assigned to the
nearest cluster):

```julia
# A conformer name from the ensemble.
cluster_conformers("/data/alphaconformers/1AKZ_A"; query = "cluster_1/.../model_1.pdb")

# A structure file; here it overrides the first-reference default.
cluster_conformers("/data/alphaconformers/1AKZ_A", "refs/1AKZ_A.pdb", "refs/1SSP_E.pdb"; query = "/path/to/my_query.pdb")
```

### Output structure

The clustering results are written into a `conformer_clustering/` subfolder of the
run folder, so they stay separate from the input `cluster_*/` predictions:

```text
/data/alphaconformers/1AKZ_A/            # the run folder you pass in
├── cluster_1/ … cluster_N/             # INPUT: AlphaConformers predictions (untouched)
└── conformer_clustering/               # OUTPUT: everything below is written here
    |
    |   # --- flat tables ---
    |-- aligned_clustering_results.csv   # Name, KMeans_K<k> (labels 1..k)
    |-- aligned_cluster_rmsd.csv         # Method, Cluster, Size, Mean_RMSD_Angstrom
    |-- agglomerative_assignments.csv    # Name, KMeans_K<k>, Sub_Cluster, Mini_Cluster (sub-clustered conformers)
    |-- dist_external_rmsds.csv          # Name, RMSD_<label>…            (only with >=1 reference)
    |-- reference_clusters.csv           # Reference, Path, Cluster, Mean_RMSD_Angstrom  (>=1 ref)
    |
    |   # --- figures ---
    |-- cluster_overview.png             # PCA scatter, colored by KMeans cluster
    |-- reference_scatter.png            # adaptive: ref-vs-ref (2) / embedding-vs-ref (1) / embedding (0)
    |-- query_path.png                   # when a reference is given (tracks the 1st) or a `query` is passed
    |
    |   # --- one folder per occupied (or surviving) KMeans cluster C ---
    |-- kmeans<C>/
    |   |-- output_cluster_members.csv   # Name, KMeans_K<k>, Sub_Cluster, Mini_Cluster
    |   |-- dendrogram.png
    |   |-- subcluster_scatter.png
    |   `-- minicluster_graph.png        # skipped if the cluster has a single sub-cluster
    |
    `-- deepaccnet/                      # only when a score table is supplied
        |-- surviving/surviving.csv      # Name, KMeans_K<k>  (the configured 25%/75% cut)
        |-- rmsd_scatter_by_score.png
        `-- kmeans_deepaccnet_filtered.png
```

Stage notes:

- Stage 1 (KMeans) and Stage 3 (agglomerative sub-clustering) always run, so the
  flat tables and the `kmeans<C>/` folders are always produced.
- The `dist_external_rmsds.csv` / `reference_clusters.csv` tables and the
  reference-aware figures appear only when at least one reference is passed.
- `deepaccnet/` and the surviving-cluster table appear only when a score table
  is supplied. With the filter on, only surviving clusters become `kmeans<C>/`
  folders, so a run with DeepAccNet typically has fewer cluster folders than the
  same run without it.

The matching fields on the returned named tuple are `results_dir` (the
`conformer_clustering/` folder above), `cluster_dirs` (the `kmeans<C>/` paths),
`deepaccnet_dir`, and `figures` (the list of written figure paths).

## Container Setup

### ColabFold Container Setup

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

### DeepAccNet Container Setup

The DeepAccNet container is an external pipeline dependency. Download the
required `deepaccnet.sif` file from Zenodo:

https://zenodo.org/records/21192188

Change `CONTAINER_DIR` to the folder where you want to store the image.

```julia
using Downloads

CONTAINER_DIR = "/path/to/container/files"
DEEPACCNET_SIF = joinpath(CONTAINER_DIR, "deepaccnet.sif")

mkpath(CONTAINER_DIR)

Downloads.download(
    "https://zenodo.org/records/21192188/files/deepaccnet.sif?download=1",
    DEEPACCNET_SIF,
)
```

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

## Citation

If you use AlphaConformers, please cite:

```bibtex
@article{Daniel2026.08.18.745512,
  author       = {Daniel, Julie and Vitoriano De Queiroz Lira, Lucas and Zea, Diego Javier},
  title        = {AlphaConformers: Structure-guided sampling enables prediction of multiple protein conformations},
  elocation-id = {2026.08.18.745512},
  year         = {2026},
  doi          = {10.64898/2026.08.18.745512},
  publisher    = {Cold Spring Harbor Laboratory},
  URL          = {https://www.biorxiv.org/content/early/2026/08/18/2026.08.18.745512},
  eprint       = {https://www.biorxiv.org/content/early/2026/08/18/2026.08.18.745512.full.pdf},
  journal      = {bioRxiv}
}
```
