# Implementation of the Hobohm algorithm I for clustering protein structures using
# the RMSD of the Cα atoms measured by USalign. We use USalign as it is slightly
# faster than the TMalign implementation.

# The sorting step will not be performed, but we will use the query structure as the 
# first cluster seed.

"""
Arguments for USalign:

 - `-mol prot`: Align proteins only 
 - `-TMcut 0.5`: Skip alignment if TM-score is unlikely to reach 0.5 
 - `-fast`: Use faster but slightly less accurate alignment method 
 - `-outfmt 2`: Use tabular format for compact output that can be parsed with CSV.jl 
"""
const _USalign_ARGS = `-mol prot -TMcut 0.5 -fast -outfmt 2`

"Parse the table output format of USalign (`-outfmt 2`) into a DataFrame"
function _read_usalign_output_table(file)
    df = DataFrames.DataFrame(CSV.File(file, delim = '\t'))
    DataFrames.rename!(df, "#PDBchain1" => :PDBchain1)
end

function _save_list(file, pdb_list, pdb_folder)
    open(file, "w") do io
        for pdb in pdb_list
            file = joinpath(pdb_folder, pdb)
            if isfile(file)
                println(io, pdb)
            else
                @warn "$file does not exist, it wont be included in the comparison."
            end
        end
    end
end

"""
    usalign(pdb_file_a, pdb_file_b)
    usalign(pdb_file_a, pdb_folder, pdb_list)

Performs structural alignment of protein structures using the USalign algorithm. The
alignment is performed using the `-fast` option, which is slightly less accurate but 
faster than the default alignment method. The `-TMcut 0.5` option is used to skip 
alignments that are unlikely to reach a TM-score of 0.5.

## Function signatures

- `usalign(pdb_file_a, pdb_file_b)` calculates the structural alignment between two 
  protein structures. It takes the path to the two protein structures in PDB format 
  as inputs.

- `usalign(pdb_file_a, pdb_folder, pdb_list)` calculates the structural alignment 
  between one protein structure and multiple other protein structures. It takes the 
  path to the protein structure in PDB format to be aligned (`pdb_file_a`), the path 
  to the directory containing the protein structures to align with (`pdb_folder`), 
  and a list of protein structures (in PDB format) contained in `pdb_folder` 
  (`pdb_list`) as inputs. You can use the `create_pdb_folder` function to create the
  `pdb_folder` using a list of target protein structures.

In both cases, the function asserts the existence of the input file(s) before proceeding 
with the alignment and returns a DataFrame with the results of the alignment.

## Implementation

This function uses the USalign tool, interfaced via the `USalign_jll` package. The actual 
alignment is performed by running the tool in a subprocess, and the output is parsed 
from the tool's tabular output format (`-outfmt 2`).

"""
function usalign(pdb_file_a, pdb_file_b)
    @assert isfile(pdb_file_a) "File $pdb_file_a does not exist"
    @assert isfile(pdb_file_b) "File $pdb_file_b does not exist"
    mktemp() do tmp_path, _
        run(pipeline(
            `$(USalign_jll.USalign()) $_USalign_ARGS $pdb_file_a $pdb_file_b`, 
            stdout=tmp_path))
        _read_usalign_output_table(tmp_path)
    end
end

function usalign(pdb_file_a, pdb_folder, pdb_list)
    @assert isfile(pdb_file_a) "File $pdb_file_a does not exist"
    @assert isdir(pdb_folder) "Folder $pdb_folder does not exist"
    # pdb_folder should end with a slash
    if !endswith(pdb_folder, '/')
        pdb_folder *= '/'
    end
    mktemp() do list, _
        _save_list(list, pdb_list, pdb_folder)
        mktemp() do output, _
            run(pipeline(
                `$(USalign_jll.USalign()) $_USalign_ARGS $pdb_file_a -dir2 $pdb_folder $list`, 
                stdout=output))
            _read_usalign_output_table(output)
        end
    end
end
