#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=test_set_development-%j.out
=#
import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

# Load necessary packages on all workers
using Distributed
using DataFrames, CSV
using Dates
using MIToS.PDB
import MIToS
using AlphaConformers
using BioStructures
using DataStructures
using Statistics

################################################################ Functions ####################################################

# Function to find the apo and holo structures for a given Uniprot ID
"""
It reads the PDB files for apo and holo structures, performs structural alignment,
and returns the apo and holo structures with the highest RMSD.
It returns the best RMSD, apo PDB ID, apo chain, holo PDB ID, and holo chain.
If no apo/holo pair is found, it returns nothing.
"""
function found_apo_holo(df_uniprot,pdb_folder)

    #Get the Apo and Holo structures
    with_ligands = filter(row -> !ismissing(row.LIGANDS), df_uniprot)
    without_ligands = filter(row -> ismissing(row.LIGANDS), df_uniprot)

    best_RMSD=nothing
    apo_to_keep=nothing
    apo_chain_to_keep=nothing
    holo_to_keep=nothing
    holo_chain_to_keep=nothing

    #For each holo structure
    for holo in eachrow(with_ligands)
        @info "Holo structure: $(holo.PDB) $(holo.CHAIN)"
        # Read the PDB file for holo structure
        pdbfilepath_holo=joinpath(pdb_folder,uppercase(String(holo.PDB))*".pdb" )
        residues_1ivo_holo = read_file(pdbfilepath_holo, PDBFile,chain=String(holo.CHAIN))

        #For each apo structure
        for apo in eachrow(without_ligands)
            @info "Apo structure: $(apo.PDB) $(apo.CHAIN)"
            # Read the PDB file for apo structure
            pdbfilepath_apo=joinpath(pdb_folder,uppercase(String(apo.PDB))*".pdb" )
            residues_1ivo_apo = read_file(pdbfilepath_apo, PDBFile,chain=String(apo.CHAIN))
            
            #get the alignment 
            try
                _,_,_,rmsd,_,_=AlphaConformers.structural_alignment(residues_1ivo_holo, residues_1ivo_apo)
                if isnothing(best_RMSD) || rmsd > best_RMSD # We want the largest RMSD, so we keep the highest one
                    best_RMSD = rmsd
                    apo_to_keep = String(apo.PDB)
                    apo_chain_to_keep = String(apo.CHAIN)
                    holo_to_keep = String(holo.PDB)
                    holo_chain_to_keep = String(holo.CHAIN)
                end
            catch e
                @warn "Error in the structural alignment: $e"
                continue
            end
            
        end
    end
    return best_RMSD, apo_to_keep, apo_chain_to_keep, holo_to_keep, holo_chain_to_keep
end

############################################################ MAIN ##############################################################
"""
This script develops a test set of apo and holo structures for AlphaConformers.
It reads a CSV file containing PDB information, processes the data to find apo and holo structures,
and saves the results in a specified output directory.
Input:
- pdb_folder: Path to the folder containing PDB files.
- output_dir: Path to the output directory where the test set will be saved.
- path_input: Path to the input CSV file containing PDB information.
Output:
- A CSV file named "info_bdd_set.csv" containing detailed information about the test set,
  including PDB IDs, chains, and RMSD values.
"""
########################## Information to fill #################################
# Define the PDB folder path
pdb_folder= abspath("/alpha/database/pdb", "pdb_files")

#Define the output directory
output_dir = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/bdd_set"

# Define the path to the input CSV file
path_input="pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek_final.csv" #File output by create_databases.jl script
################################################################################

@info "Starting test set development script"
@show "Date de debut ", Dates.format(now(), "HH:MM:SS")

#Create the output directory if it does not exist
if !isdir(output_dir)
    mkdir(output_dir)
end

# Read the input CSV file into a DataFrame
df_input=DataFrames.DataFrame(CSV.File(path_input,
comment="#", missingstring=["", "None"]))

# Initialize the DataFrame to store the test set information
info_test_set = DataFrame(PDB_apo=String[], CHAIN_apo=String[], PDB_holo=String[], CHAIN_holo=String[],RMSD=Float64[])

# Get unique Uniprot IDs
unique_uniprot = unique(df_input.UNIPROT)
@show "Number of unique uniprot: ", length(unique_uniprot)
# Process each Uniprot ID
for uniprot in unique_uniprot
    @info "Processing uniprot : $uniprot"
    df_uniprot = filter(row -> row.UNIPROT == uniprot, df_input)
    
    #get the apo and holo structures
    result = found_apo_holo(df_uniprot, pdb_folder)

    # Check if the result is valid
    # If result is nothing or contains any nothing values, skip to the next Uniprot
    if result === nothing || any(isnothing, result)
        @warn "No apo and holo structure found for $uniprot"
        continue
    end
    best_RMSD, apo_to_keep, apo_chain_to_keep, holo_to_keep, holo_chain_to_keep = result
    @info "Best RMSD for $uniprot: $best_RMSD"
    @show apo_to_keep
    @show holo_to_keep

    #Save the best apo and holo structures
    pdbfilepath_apo=joinpath(pdb_folder,uppercase(apo_to_keep)*".pdb" )
    pdb_file = read_file(pdbfilepath_apo, PDBFile,chain=apo_chain_to_keep) # Read the PDB file for apo structure
    pdb_file_path = joinpath(output_dir,"$(apo_to_keep)_$(apo_chain_to_keep).pdb") 
    write_file(pdb_file_path, pdb_file, MIToS.PDB.PDBFile) # Save the apo structure to the output directory

    pdbfilepath_holo=joinpath(pdb_folder,uppercase(holo_to_keep)*".pdb" )
    pdb_file = read_file(pdbfilepath_holo, PDBFile,chain=holo_chain_to_keep) # Read the PDB file for holo structure
    pdb_file_path = joinpath(output_dir,"$(holo_to_keep)_$(holo_chain_to_keep).pdb")
    write_file(pdb_file_path, pdb_file, MIToS.PDB.PDBFile) # Save the holo structure to the output directory
    
    # Add the information to the DataFrame
    @info "Ajout des informations au DataFrame"
    push!(info_test_set, (PDB_apo=apo_to_keep, CHAIN_apo=apo_chain_to_keep, PDB_holo=holo_to_keep, CHAIN_holo=holo_chain_to_keep, RMSD=best_RMSD))
    
end
# Save the DataFrame to a CSV file
info_test_set_path= joinpath(output_dir, "info_bdd_set.csv")
CSV.write(info_test_set_path, info_test_set)

@info "End"
############################################################ End ##############################################################