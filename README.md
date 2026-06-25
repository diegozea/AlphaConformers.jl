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

Run the default pipeline with input preparation, ColabFold prediction, and output
triage:

```julia
using AlphaConformers

query_struct = "1ABC_A.pdb"
pdb_folder = "datasets/pdb/mmcif_files"
output_dir = "outputs/1ABC_A"

foldseek_dbs = [
    "datasets/foldseek/fullpdb",
    "datasets/foldseek/afdb",
]

sif_path = "containers/ColabFold_AF2_1-5-5.sif"
cache_dir = "cache/colabfold"

alphaconformers(
    sif_path,
    cache_dir;
    query_struct,
    pdb_folder,
    output_dir,
    databases = foldseek_dbs,
    n_threads = 16,
    evalue_cutoff = 1e-5,
    cutoff = 1.0,
)
```

The query file name should include the PDB code and chain, for example
`1ABC_A.pdb`.

The Foldseek database paths in `databases` must be set for your machine. AlphaConformers
does not assume a default local database path. When several databases are used,
exactly one database name should contain `pdb`, such as `fullpdb`; that lets
AlphaConformers pass the PDB Foldseek result folder to triage instead of using
database order.

## Pipeline Steps

Each step can also be run separately with `alphaconformers`. The examples below
reuse `query_struct`, `pdb_folder`, `output_dir`, and `foldseek_dbs` from the
quick start.

```julia
# Prepare only. Requires query_struct, pdb_folder, output_dir, and databases.
alphaconformers(;
    query_struct,
    pdb_folder,
    output_dir,
    databases = foldseek_dbs,
    predict = false,
    triage = false,
)
```

Prediction uses the default `structure_predictor = run_alphafold`. That predictor
expects two positional arguments before the semicolon: `sif_path` and `cache_dir`.
`sif_path` is the ColabFold Apptainer `.sif` image:

```julia
sif_path = "containers/ColabFold_AF2_1-5-5.sif"
```

The `cache_dir` folder is important for ColabFold runs. It is mounted inside the
Apptainer container as `/cache/colabfold`, where ColabFold looks for its model
weights and cache files. Choose a permanent location and reuse the same folder
across runs so the weights do not need to be downloaded again. This saves time
and avoids repeated large downloads. The cache can require around 4 GB of disk
space or more. Prediction outputs are not written to this folder; they are
written inside each cluster's `af/` directory.

```julia
cache_dir = "cache/colabfold"
```

Now prediction can be run by itself:

```julia
# Predict only. Requires output_dir plus the predictor-specific arguments.
alphaconformers(
    sif_path,
    cache_dir;
    output_dir,
    prepare = false,
    triage = false,
)
```

After prediction, triage can also be run separately. It requires `output_dir`,
`query_struct`, and a SIFTS mapping:

```julia
sifts_uniprot_mapping = get_uniprot_mapping()

# Triage only. Requires output_dir, query_struct, and a SIFTS mapping.
alphaconformers(;
    output_dir,
    query_struct,
    prepare = false,
    predict = false,
    triage = true,
    sifts_uniprot_mapping,
)
```

When preparation and triage run in the same `alphaconformers` call, the
Foldseek result folder is passed between steps automatically. If preparation uses
several databases, AlphaConformers chooses the unique database whose name contains
`pdb`, ignoring case. If it cannot choose one, it stops before running Foldseek
and asks for `foldseek_results_folder`. For resumed runs, triage searches
`output_dir` for a single `*_results` folder containing `.m8` files. If several
such folders exist, pass the folder explicitly:

```julia
alphaconformers(;
    output_dir,
    query_struct,
    prepare = false,
    predict = false,
    triage = true,
    sifts_uniprot_mapping,
    foldseek_results_folder = joinpath(output_dir, "fullpdb_results"),
)
```

Any extra positional arguments and any unknown keyword arguments are passed only
to `structure_predictor`. If `predict=false`, predictor-specific arguments are
rejected so they are not silently ignored.

For a custom predictor:

```julia
alphaconformers(
    predictor_arg;
    output_dir,
    prepare = false,
    triage = false,
    structure_predictor = my_predictor,
    predictor_keyword = value,
)
```

The lower-level step functions remain available:

```julia
prepared = prepare_inputs(query_struct, pdb_folder, output_dir; databases = foldseek_dbs)
run_alphafold(output_dir, sif_path, cache_dir)
triage_outputs(
    output_dir,
    query_struct,
    sifts_uniprot_mapping;
    foldseek_results_folder = prepared.foldseek_results_folder,
)
```

`prepare_inputs` returns a small `PreparedInputs` object with the query path, the
output folder, and the Foldseek result folder. You can ignore this return value
when you only need the files written to `output_dir`.

`triage_outputs` is the final pipeline step. It currently filters predictions
through `found_best_prediction`; future versions may add more output triage on
top of that step.

AlphaConformers also provides helpers for AlphaFold3 and Boltz2:

```julia
run_alphafold3(output_dir, sif_image_dir, model_parameters_dir, db_dir)
run_boltz2(output_dir, boltz_sif_path, boltz_cache_dir)
```

These helpers expect local container images and local model or database paths.
They do not download model weights.

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
separate from `cache_dir`, which stores reusable ColabFold files and can require 
around 4 GB or more.


## Conformer Clustering

`cluster_conformers` groups an AlphaConformer prediction ensemble for one
protein system by shape. It runs three stages - KMeans clustering, an optional
DeepAccNet score filter, and agglomerative sub-clustering, it writes a
per-system output tree plus the essential figures. 

```julia
using AlphaConformers

# Reference-free: just cluster the ensemble under data_root/<system>.
cluster_conformers("1AKZ")

# One reference (apo) or two references (apo + holo). The reference mode is
# auto-detected from the stems; there is no mode flag and no chain argument.
cluster_conformers("1AKZ", "1AKZ_A")
cluster_conformers("1AKZ", "1AKZ_A", "1SSP_E")

# Common options (shown with their defaults).
cluster_conformers("1AKZ", "1AKZ_A", "1SSP_E";
    data_root = "/data/alphaconformers/pdb",  # holds <system>/ AlphaConformers output trees
    refs_dir  = "refs",                        # holds the apo/holo reference files
    out_dir   = "results",                     # output root
    kmeans_k  = 30,
    threshold = 1.0,                           # agglomerative cut in Å
    linkage   = :average,
    seed      = 42)
```

The input under `data_root/<system>/` is an AlphaConformers output tree. Only the
predicted structures are clustered: the `.pdb`/`.cif` files inside `models/`
folders (the prediction location, `cluster_*/.../predictions/sequences/models/`).

A reference stem is `<PDBID>_<CHAIN>` (for example `1AKZ_A`) resolved against
`refs_dir`.

To run the optional DeepAccNet score filter (Stage 2), pass a score table - a
`DataFrame` or a path to a CSV with a `Name` column and the score column:

```julia
cluster_conformers("1AKZ", "1AKZ_A", "1SSP_E"; score_table = "deepaccnet_results.csv")
```

### Query path

The query-path figure (`query_path.png`) ranks the clusters by their RMSD to a
**query** structure, starting from the cluster the query falls into. The query
defaults to a reference, so you usually do not set it:

- One reference (apo only): that reference is used as the query.
- Two references (apo + holo): the **apo** reference is used as the query.
- Reference-free: no query is used unless you pass one, so no query-path figure.

```julia
# query_path.png is produced automatically - the apo reference is the query.
cluster_conformers("1AKZ", "1AKZ_A")            # 1 reference  → that ref is the query
cluster_conformers("1AKZ", "1AKZ_A", "1SSP_E")  # 2 references → apo is the query
```

To override the default, pass `query` explicitly - either a conformer name (as it
appears in the clustering table) or a path to a structure file (assigned to the
nearest cluster):

```julia
# A conformer name from the ensemble.
cluster_conformers("1AKZ"; query = "cluster_1/.../model_1.pdb")

# A structure file; here it overrides the apo default in two-reference mode.
cluster_conformers("1AKZ", "1AKZ_A", "1SSP_E"; query = "/path/to/my_query.pdb")
```

### Output structure

Everything for one run is written under `out_dir/<system>/`. The tree has three
parts: flat tables and figures at the system root, one `kmeans<C>/` folder per
cluster, and an optional `deepaccnet/` folder when a score table was supplied.

```text
results/<SYSTEM>/
|
|   # --- system-root tables ---
|-- aligned_clustering_results.csv   # Stage 1 KMeans labels.  Columns: Name, KMeans_K<k> (labels 0..k-1)
|-- aligned_cluster_rmsd.csv         # Stage 1 intra-cluster RMSD.  Columns: Method, Cluster, Size, Mean_RMSD_Angstrom
|-- agglomerative_assignments.csv    # Stage 3 flat table for ALL conformers.  Columns: Name, KMeans_K<k>, Sub_Cluster, Mini_Cluster
|-- dist_external_rmsds.csv          # conformer-to-reference RMSD (only with >=1 reference).  Columns: Name, RMSD_Apo[, RMSD_Holo]
|-- reference_clusters.csv           # each reference's best cluster (only with >=1 reference).  Columns: Reference, Stem, Cluster, Mean_RMSD_Angstrom
|
|   # --- system-level figures ---
|-- cluster_overview.png             # PCA scatter, colored by KMeans cluster
|-- reference_scatter.png            # adaptive: apo-vs-holo (2 refs) / embedding-vs-reference (1 ref) / plain embedding (0 refs)
|-- query_path.png                   # when a reference is given (tracks apo) or a `query` is passed
|
|   # --- one folder per occupied KMeans cluster C ---
|   # (when the score filter ran, only the SURVIVING clusters get a folder)
|-- kmeans<C>/
|   |-- members.csv                  # this cluster's members.  Columns: Name, KMeans_K<k>, Sub_Cluster, Mini_Cluster
|   |-- dendrogram.png               # agglomerative tree of this cluster's members
|   |-- subcluster_scatter.png       # members colored by sub-cluster
|   `-- minicluster_graph.png        # sub-cluster MST graph (skipped if the cluster has one sub-cluster)
|
`-- deepaccnet/                      # Stage 2, only when a score table is supplied
    |-- surviving/
    |   `-- surviving_rel<rel>_frac<frac>.csv   # the configured surviving cut (default rel=25, frac=75).  Columns: Name, KMeans_K<k>
    |-- rmsd_scatter_by_score.png    # conformer RMSD colored by DeepAccNet score
    `-- kmeans_deepaccnet_filtered.png   # clusters with the filtered-out ones de-emphasized
```

Stage notes:

- Stage 1 (KMeans) and Stage 3 (agglomerative sub-clustering) always run, so the
  three system-root tables and the `kmeans<C>/` folders are always produced.
- The `dist_external_rmsds.csv` / `reference_clusters.csv` tables and the
  reference-aware figures appear only when at least one reference is passed.
- `deepaccnet/` and the surviving-cluster tables appear only when a score table
  is supplied. With the filter on, only surviving clusters become `kmeans<C>/`
  folders, so a run with DeepAccNet typically has fewer cluster folders than the
  same run without it.

The matching fields on the returned named tuple are `system_dir` (the root
above), `cluster_dirs` (the `kmeans<C>/` paths), `deepaccnet_dir`, and `figures`
(the list of written figure paths).


## Common Helpers

Some lower-level helpers are useful when building custom workflows:

- `run_foldseek`: run Foldseek and save hit tables, alignments, and aligned
  structures.
- `merge_tables`: merge Foldseek hit tables.
- `merge_msas`: merge Foldseek alignments.
- `get_uniprot_mapping`: read or download SIFTS UniProt mappings.
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
