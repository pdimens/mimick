"""
The main function that takes arguments and does the thing
"""
function mimick(prefix::String, fasta_files::Vector{String}, coverage::Float64, mol_cov::Float64, mol_len::Int64, read_distance::Int, distance_stdev::Int, read_len::Vector{Int}, singletons::Float64, circular::Bool, attempts::Int)
    schemas = setup_schema(fasta_files, coverage, read_len)
    params = SimParams(prefix, 0.001, read_distance, distance_stdev, read_len[1], read_len[2], mol_len, mol_cov, singletons, circular = circular, attempts = attempts)
    
    
    return schemas, params
end