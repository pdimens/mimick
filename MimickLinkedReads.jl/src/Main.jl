"""
The main function that takes arguments and does the thing
"""
function mimick(fasta_files::Vector{String}; format::String, prefix::String = "simulated/SIM", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = 0)
    if seed > 0
        Random.seed!(seed)
    end
    schema = setup_schema(fasta_files, coverage, read_len)
    barcodes = setup_barcodes(Symbol(format))
    params = SimParams(prefix, 0.001, insert_length, insert_stdev, read_len[1], read_len[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
    mkpath(dirname(prefix))
    writer = FastqWriter(prefix, Symbol(format))

    while !isempty(schema)
        candidates = keys(schema)
        n_mol = get_n_molecules(params)
        molecule_targets = rand(candidates, n_mol)
        for target in molecule_targets
            molsize = get_molecule_size(params, length(schema[target].sequence))
            frags = calculate_fragments(params, molsize)
            molecule = get_sequences(schema[target], params, get_next!(barcodes), molsize, frags)
            Base.Threads.atomic_add!(schema[target].tracker.reads_current, length(frags))
            submit!(writer, molecule)
        end
        update!(schema)
    end
    stop!(writer)
end