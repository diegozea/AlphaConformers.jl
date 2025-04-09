#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=4
#PBS -l mem=16g
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

using DataFrames, CSV
using Dates
using MIToS.PDB
using AlphaConformers
using MIToS.PDB


function foldseek_similar_pdb(df_final,FOLDSEEK_DB,pdb_folder)
    database=DataFrame(UNIPROT=String[],PDB=String[],NB_SIMILAR_PROT=Int64[]) # Create an empty DF
    df_merged=DataFrame()
    for row in eachrow(df_final)
        uniprot=row.UNIPROT # Check for every uniprot accession 
        println(uniprot)
        #Check that we haven't already done this
        is_present= uniprot in database.UNIPROT
        if !is_present
            df_uniprot = filter(row -> row.UNIPROT == uniprot, df_final)# To focus on the uniprot 
            #Look if we have both conformation for this uniprot
            is_apo = is_apo = any(ismissing, df_uniprot.LIGANDS) && any(!ismissing, df_uniprot.LIGANDS)
            if is_apo 
                df_uniprot_res = filter(row -> ismissing(row.RESOLUTION) || row.RESOLUTION < 2.5, df_uniprot) #Pour recuperer fichier avec bonne résolution
                df_uniprot_res_date = filter(row -> row.DATE_RELEASE >= Date(2018-04-30), df_uniprot_res) #get pdb release after the training of AF2
                if !isempty(df_uniprot_res_date) #check that we still have row
                    is_apo = is_apo = any(ismissing, df_uniprot_res_date.LIGANDS) && any(!ismissing, df_uniprot_res_date.LIGANDS)
                    if is_apo #Check that we still have apo and holo conformation 

                        #Check that the apo and holo conformation have at least 1A of differences 
                        holo = filter(row -> !ismissing(row.LIGANDS), df_uniprot_res_date)
                        apo = filter(row -> ismissing(row.LIGANDS), df_uniprot_res_date)
                        for row_holo in eachrow(holo)
                            pdb_holo=row_holo.PDB
                            pdbfilepath_holo=joinpath(pdb_folder,uppercase(pdb_holo)*".pdb" ) # Use the temporary file
                            if isfile(pdbfilepath_holo) 
                                #Get the right chain only
                                residues_1ivo_holo = read(pdbfilepath_holo, PDBFile)
                                chain_residues_pdb_holo = MIToS.PDB.residuesdict(residues_1ivo_holo, "1", String(row_holo.CHAIN), "ATOM") #Returns dictionary with position => residues details 
                                rmsd_list = []
                                for row_apo in eachrow(apo)
                                    pdb_apo=row_apo.PDB
                                    pdbfilepath_apo=joinpath(pdb_folder,uppercase(pdb_apo)*".pdb" ) # Use the temporary file
                                    if isfile(pdbfilepath_apo) 
                                        #Get the right chain only
                                        residues_1ivo_apo = read(pdbfilepath_apo, PDBFile)
                                        chain_residues_pdb_apo = MIToS.PDB.residuesdict(residues_1ivo_apo, "1", String(row_apo.CHAIN), "ATOM") #Returns dictionary with position => residues details 
                                        
                                        #get the alignment 
                                        common_keys = Set(keys(chain_residues_pdb_holo)) ∩ Set(keys(chain_residues_pdb_apo))
                                        positions = SortedSet(parse.(Int,common_keys))  # Retrieves common positions (same key in all files)
                                        # Extracting residues for each file as an ordered list
                                        residues_holo = [chain_residues_pdb_holo[string(pos)] for pos in positions] 
                                        residues_apo = [chain_residues_pdb_apo[string(pos)] for pos in positions]
                                        try  
                                            # Superposition and calculation of RMSD
                                            _, _,rmsd = superimpose(residues_holo, residues_apo)
                                        catch e 
                                            @show common_keys
                                            @warn "❌ Error while superimposing structures "
                                            rmsd = missing # we will not annalyse those pdb 
                                        end 
                                        push!(rmsd_list, rmsd)
                                    else 
                                        @warn "The PDB file of $pdb_apo doesn't exist"
                                        index_to_remove = findall(row -> row.PDB == pdb_apo, df_uniprot_res_date)
                                        delete!(df_uniprot_res_date, index_to_remove)
                                    end
                                end
                                if mean(rmsd_list)<1
                                    index_to_remove = findall(row -> row.PDB == pdb_holo && row.CHAIN == String(row_holo.CHAIN), df_uniprot_res_date)
                                    delete!(df_uniprot_res_date, index_to_remove)
                                end
                            else 
                                @warn "The PDB file of $pdb_holo doesn't exist"
                                index_to_remove = findall(row -> row.PDB == pdb_holo, df_uniprot_res_date)
                                delete!(df_uniprot_res_date, index_to_remove)
                            end 
                        end 
                        is_apo = is_apo = any(ismissing, df_uniprot_res_date.LIGANDS) && any(!ismissing, df_uniprot_res_date.LIGANDS)
                        if is_apo 
                            for row in eachrow(df_uniprot_res_date)
                                pdb_file=row.PDB
                                pdbfilepath=joinpath(pdb_folder,uppercase(pdb_file)*".pdb" ) # Use the temporary file
                                if isfile(pdbfilepath) #Check that the file is not already downloaded
                                    out_file=AlphaConformers.foldseek_search(pdbfilepath,db_path=FOLDSEEK_DB)
                                    df_result=AlphaConformers.read_foldseek_search_results(out_file)
                                    #Filter the evalue to keep the one < 1e-3
                                    df_result_evalue = filter(row.evalue < 1e-3, df_result)
                                    push!(database,(uniprot,pdb_file,size(df_result_evalue)[1]))
                                else 
                                    @warn "The PDB file of $pdb_file doesn't exist"
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        df_merged = innerjoin(database, df_final, on=[:UNIPROT, :PDB])
    end
    return df_merged
end

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
pdb_folder= abspath("/alpha/database/pdb", "pdb_files")

#Get the file with the result 
file_path_df_final="pdb_information_details_final_1000.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_final,20))
println(size(df_final))

df_merged=foldseek_similar_pdb(df_final,FOLDSEEK_DB,pdb_folder)
println(first(df_merged,20))
println(size(df_merged))