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