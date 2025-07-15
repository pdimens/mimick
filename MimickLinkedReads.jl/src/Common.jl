"""
`safe_read(filename::String)`

Read in a file that might be gzipped. Returns an IO stream.
"""
function safe_read(filename::String)
    try
        readline(GzipDecompressorStream(open(filename, "r")))
        return GzipDecompressorStream(open(filename, "r"))
    catch e
        if isa(e, CodecZlib.ZlibError) || isa(e, EOFError)
            return open(filename, "r")
        else
            rethrow(e)
        end
    end
end

"""
`is_complete(schema::Schema)`

Returns `true` if the current number of reads is equal to or greater than the
number of reads required.
"""
is_complete(schema::Schema)::Bool = schema.tracker.reads_current >= schema.tracker.reads_required

"""
`update!(schemas::Dict{String, Schema})`

Parses through the Dict of Schema and deletes entries where
the number of generated reads is greater than or equal to the
reads required for that contig/interval.
"""
function update!(schemas::Dict{String, Schema})
    for i in keys(schemas)
        if is_complete(schemas[i])
            delete!(schemas, i)
        end
    end
end

"""
`get_n_molecules(::SimParams)`

Randomly sample a `Distribution` (or `UnitRange`) and return a rounded `Int`.
"""
@inline get_n_molecules(params::SimParams)::Int = get_n_molecules(params.n_molecules)
@inline get_n_molecules(dist::Distribution{Univariate, Continuous})::Int = trunc(Int,round(rand(dist), digits = 0))
@inline get_n_molecules(dist::UnitRange{Int})::Int = dist.start