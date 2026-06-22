```@meta
CurrentModule = AlphaConformers
```

# Database Setup

AlphaConformers searches local protein structure databases with Foldseek. The database is
built once from structure files, then reused for many searches.

You do not need to install Foldseek separately when you use AlphaConformers from Julia.
The package uses `Foldseek_jll`, which provides the Foldseek executable.

## Use the Full PDB

For AlphaConformers, use the full PDB archive rather than a reduced or representative
set such as PDB100. Different PDB entries for the same or related proteins can capture
different conformational states. A representative database can remove those alternatives.

Download the full PDB in mmCIF format when possible. The wwPDB and RCSB document the
official archive download locations and recommend `rsync` for maintaining local copies:

- [wwPDB PDB Archive Downloads][wwpdb-downloads]
- [RCSB PDB File Download Services][rcsb-downloads]

Foldseek can read common structure formats such as `.pdb`, `.cif`, `.pdb.gz`, and
`.cif.gz`. The Foldseek manual has more details:

- [Foldseek documentation][foldseek-docs]

## Prepare the Structure Folder

Put the structures in one folder:

```text
datasets/pdb/mmcif_files/
    100D.cif
    101M.cif
    ...
```

Use storage with enough free space. For a shared cluster, choose a path that will still
be available when other jobs run later.

## Create the Foldseek Database

Start Julia from the AlphaConformers project:

```bash
julia --project=.
```

In the Julia REPL, load Foldseek and define the paths once:

```julia
import Foldseek_jll

STRUCTURE_PATH = "datasets/pdb/mmcif_files"
DB_PATH = "datasets/pdb/fullpdb_mmcif_files"
TMP_PATH = "datasets/pdb/tmp"
```

`STRUCTURE_PATH` is the folder with the `.cif` or `.pdb` files. `DB_PATH` is the
Foldseek database name. It is not a folder. Foldseek will create several files that
begin with this name, for example:

```text
datasets/pdb/fullpdb_mmcif_files
datasets/pdb/fullpdb_mmcif_files.dbtype
datasets/pdb/fullpdb_mmcif_files.index
datasets/pdb/fullpdb_mmcif_files.lookup
```

Create the database:

```julia
run(`$(Foldseek_jll.foldseek()) createdb $STRUCTURE_PATH $DB_PATH`)
```

For large databases, keep the command output or log file. Foldseek reports how many
entries were skipped because they were too short or were not proteins.

## Create the Index

The index makes repeated searches faster. Create a temporary folder and build the index:

```julia
mkpath(TMP_PATH)
run(`$(Foldseek_jll.foldseek()) createindex $DB_PATH $TMP_PATH`)
```

Keep the database files together. Move or copy all files that start with `DB_PATH`, not
only the file named exactly like `DB_PATH`.

## Check That It Works

Search one protein structure against the new database:

```julia
QUERY_PATH = "datasets/pdb/mmcif_files/101M.cif"

mktempdir() do work
    output = joinpath(work, "test.m8")
    tmp = joinpath(work, "tmp")

    run(`$(Foldseek_jll.foldseek()) easy-search $QUERY_PATH $DB_PATH $output $tmp`)
    println(read(output, String))
end
```

If the command prints one or more hits, the database is usable.

## Use the Database with AlphaConformers

Pass the same `DB_PATH` to `prepare_inputs`:

```julia
using AlphaConformers

input_pdb = "1ABC_A.pdb"
pdb_folder = "datasets/pdb/mmcif_files"
output_dir = "outputs/1ABC_A"

prepare_inputs(
    input_pdb,
    pdb_folder,
    output_dir;
    db = [DB_PATH],
)
```

The query file name should include the PDB code and chain, for example `1ABC_A.pdb`.

## Updating a Database

When you add, remove, or replace many structure files, rebuild the database and recreate
the index. For a production database, write the new database to a new `DB_PATH` first,
test it, and then switch your workflows to the new path.

[foldseek-docs]: https://github.com/steineggerlab/foldseek
[rcsb-downloads]: https://www.rcsb.org/docs/programmatic-access/file-download-services
[wwpdb-downloads]: https://www.wwpdb.org/ftp/pdb-ftp-sites
