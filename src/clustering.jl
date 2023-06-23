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

function _save_list(file, pdb_list)
    open(file, "w") do io
        for pdb in pdb_list
            println(io, pdb)
        end
    end
end

function usalign_one2one(pdb_file_a, pdb_file_b)
    @assert isfile(pdb_file_a) "File $pdb_file_a does not exist"
    @assert isfile(pdb_file_b) "File $pdb_file_b does not exist"
    mktemp() do tmp_path, _
        run(pipeline(
            `$(USalign_jll.USalign()) $_USalign_ARGS $pdb_file_a $pdb_file_b`, 
            stdout=tmp_path))
        _read_usalign_output_table(tmp_path)
    end
end

function usalign_one2many(pdb_file_a, pdb_folder, pdb_list)
    @assert isfile(pdb_file_a) "File $pdb_file_a does not exist"
    @assert isdir(pdb_folder) "Folder $pdb_folder does not exist"
    # pdb_folder should end with a slash
    if !endswith(pdb_folder, '/')
        pdb_folder *= '/'
    end
    mktemp() do list, _
        _save_list(list, pdb_list)
        mktemp() do output, _
            run(pipeline(
                `$(USalign_jll.USalign()) $_USalign_ARGS $pdb_file_a -dir2 $pdb_folder $list`, 
                stdout=output))
            _read_usalign_output_table(output)
        end
    end
end
