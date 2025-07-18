"""
    safe_read(filename::String)

Read in a file that might be gzipped. Returns an IO stream.
"""
function safe_read(filename::String)
    if !isfile(filename)
        error("$file does not exist.")
    end
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
    interperet_format(fmt::String) -> Tuple{Symbol,Symbol}

Given a `fmt` (format) string, interprets the barcode type and returns
a Tuple of the barcode format and output format style (e.g. `(:tellseq, :tellseq)`).
Accepts
- "haplotagging" -> `(:haplotagging, :haplotagging)`
- "stlfr" -> `(:stlfr, :stlfr)`
- "tellseq" -> `(:tellseq, :tellseq)`
- "tenx" -> `(:tenx, :tenx)`
- "standard" -> `(:tellseq, :standard)`
- "standard:haplotagging" -> `(:haplotagging, :standard)`
- "standard:stlfr" -> `(:stlfr, :standard)`
"""
function interperet_format(fmt::String)::Tuple{Symbol,Symbol}
    if occursin("standard", fmt)
        out_format = :standard
        if occursin("stlfr", fmt)
            barcode_format = :stlfr
        elseif occursin("haplotagging", fmt)
            barcode_format = :haplotagging
        else
            barcode_format = :tellseq
        end

    elseif fmt ∈ ["stlfr", "haplotagging","tellseq","tenx"]
        out_format = Symbol(fmt)
        barcode_format = Symbol(fmt)
    else
        error("Barcode format $fmt is not recognized as one of: haplotagging, stlfr, tellseq, tenx")
    end
    return (barcode_format, out_format)
end


"""
    is_incomplete(schema::Pair{String,Schema}) -> Bool

Given a "string" => Schema Pair, returns `true` if the current number of reads is equal less than the number of reads required.
"""
is_incomplete(schema::Pair{String,Schema})::Bool = schema.second.tracker.reads_current < schema.second.tracker.reads_required

"""
    get_n_molecules(::SimParams) -> Int

Randomly sample a `Distribution` (or `UnitRange`) and return a rounded `Int`.
"""
@inline get_n_molecules(params::SimParams)::Int = get_n_molecules(params.n_molecules)
@inline get_n_molecules(dist::Distribution{Univariate, Continuous})::Int = trunc(Int,round(rand(dist), digits = 0))
@inline get_n_molecules(dist::UnitRange{Int})::Int = dist.start