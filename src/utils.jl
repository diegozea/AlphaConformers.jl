"""
    _create_empty_folder(path)

Helper function to create an empty folder at the given path. If the folder already exists,
all its contents will be deleted.
"""
function _create_empty_folder(path)
    if isdir(path)
        @warn "The folder $path already exists; all its contents will be deleted."
        rm(path, recursive=true, force=true)
    end
    @info "Creating folder $path"
    mkdir(path)
end


# --- small helpers -----------------------------------------------------------

# Parse the numeric part of a PDB residue number like "34A" -> 34
_auth_num(resnum::AbstractString) = something(
    tryparse(Int, replace(resnum, r"[A-Za-z]" => "")),
    typemax(Int),
)

# Parse insertion code like "34A" -> 'A', "34" -> ' '
_ins_char(resnum::AbstractString) = begin
    m = match(r"[A-Za-z]$", resnum)
    m === nothing ? ' ' : m.match[1]
end
# --- amino acid name test (3-letter codes) -----------------------------

const AA3 = Set([
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    # common variants seen by AlphaFold
    "MSE","SEC","PYL"
])

_is_peptide_resname(resname::AbstractString) =
    uppercase(strip(resname)) in AA3

# Best-effort element guess if missing (common for CA-only traces)
function _guess_element(atomname::AbstractString)
    s = uppercase(strip(atomname))
    isempty(s) && return "?"
    # Protein atom names typically start with C/N/O/S/P/H -> use first char
    if s[1] in ('C','N','O','S','P','H')
        return string(s[1])
    end
    # Otherwise allow 2-letter element symbols for ions/metal-like names
    if length(s) ≥ 2
        two = s[1:2]
        if two in ("ZN","FE","MG","NA","CL","BR","MN","CO","NI","CU","CD","HG","PB","SR","CS","BA")
            return two
        end
    end
    return string(s[1])
end
"""
    patch_mmcif_for_alphafold(infile, outfile;
                                entry_id=nothing,
                                exptl_method="computational model",
                                force_peptide_types=false,
                                fill_missing_elements=true)

Read an mmCIF (e.g. written by MIToS from a PDB) and write a "patched" mmCIF that
AlphaFold/ColabFold's `mmcif_parsing.py` can use as a template.

Returns `outfile`.
"""
function patch_mmcif_for_alphafold(infile::AbstractString, outfile::AbstractString;
    entry_id::Union{Nothing,String}=nothing,
    exptl_method::AbstractString="computational model",
    force_peptide_types::Bool=false,
    fill_missing_elements::Bool=true,
)
    # 1) Read residues from mmCIF
    residues = MIToS.PDB.read_file(infile, MIToS.PDB.MMCIFFile; label=true)

    # Sort residues to make label_seq_id / entity_poly_seq consistent
    sort!(residues, by = r -> (tryparse(Int, r.id.model) === nothing ? 1 : parse(Int, r.id.model),
                                r.id.chain, _auth_num(r.id.number), _ins_char(r.id.number)))

    # 2) Assign PDBe_number sequentially per chain (drives label_seq_id in MIToS writer)
    chain_counter = Dict{String,Int}()

    for r in residues
        c = r.id.chain
        n = get!(chain_counter, c, 0) + 1
        chain_counter[c] = n

        # Replace residue id (PDBResidueIdentifier is immutable; PDBResidue is mutable)
        r.id = MIToS.PDB.PDBResidueIdentifier(
            string(n),      # PDBe_number => becomes label_seq_id
            r.id.number,    # keep original PDB numbering => auth_seq_id
            r.id.name, r.id.group, r.id.model, r.id.chain
        )

        # Optionally fill missing element strings so type_symbol isn't "?"
        if fill_missing_elements
            r.atoms = [
                isempty(a.element) ?
                    MIToS.PDB.PDBAtom(a.coordinates, a.atom, _guess_element(a.atom),
                            a.occupancy, a.B, a.alt_id, a.charge) :
                    a
                for a in r.atoms
            ]
        end
    end

    # 3) Build atom_site twice (label + auth) then merge so we have BOTH sets of fields
    mm_label = BioStructures.MMCIFDict(residues; label=true)
    mm_auth  = BioStructures.MMCIFDict(residues; label=false)

    for (k, v) in mm_auth
        if !haskey(mm_label, k)
            mm_label[k] = v
        end
    end
    mm = mm_label

    # 4) Add minimal non-atom_site categories AlphaFold expects
    eid = entry_id === nothing ? splitext(basename(infile))[1] : entry_id
    mm["_entry.id"] = [eid]
    mm["_exptl.method"] = [String(exptl_method)]

    # REQUIRED by AlphaFold 3 templates code
    mm["_pdbx_database_status.recvd_initial_deposition_date"] = ["2021-09-30"]
    mm["_database_PDB_rev.date_original"] = ["2021-09-30"]
    mm["_pdbx_audit_revision_history.revision_date"] = ["2021-09-30"]

    # Chains present
    chains = unique([r.id.chain for r in residues])

    # struct_asym: chain id -> entity id
    chain2entity = Dict{String,String}(c => string(i) for (i, c) in enumerate(chains))
    mm["_struct_asym.id"] = chains
    mm["_struct_asym.entity_id"] = [chain2entity[c] for c in chains]

    # entity_poly_seq: (entity_id, num, mon_id) for all residues in each chain
    ent = String[]; num = String[]; mon = String[]
    for c in chains
        seqnum = 0
        for r in residues
            r.id.chain == c || continue
            seqnum += 1
            push!(ent, chain2entity[c])
            push!(num, string(seqnum))
            push!(mon, r.id.name)  # uses your already-renamed residue names
        end
    end
    mm["_entity_poly_seq.entity_id"] = ent
    mm["_entity_poly_seq.num"] = num
    mm["_entity_poly_seq.mon_id"] = mon

    # chem_comp: unique monomers and their types (AlphaFold checks for "peptide" substring)
    monomers = unique(mon)
    sort!(monomers)

    mm["_chem_comp.id"] = monomers
    mm["_chem_comp.type"] = [
        (force_peptide_types || _is_peptide_resname(m)) ?
            "L-peptide linking" : "non-polymer"
        for m in monomers
    ]

    mm["_chem_comp.name"] = monomers  # optional, but harmless

    # 5) Write patched cif
    open(outfile, "w") do io
        BioStructures.writemmcif(io, mm)
    end

    return outfile
end

"""
    read_a3m(path) -> (ids, seqs)

Parse an A3M-formatted multiple sequence alignment file.
Input : 
- `path` : path to the `.a3m` file.
Output : 
- `ids`  : vector of sequence identifiers (header lines without the leading `>`).
- `seqs` : vector of sequences, with lowercase insertion columns removed.
"""
function read_a3m(path)
    msa = MIToS.MSA.read_file(path, MIToS.MSA.A3M)
    ids = MIToS.MSA.sequencenames(msa)
    seqs = [join(msa[i, :]) for i in 1:size(msa, 1)]
    return ids, seqs
end
