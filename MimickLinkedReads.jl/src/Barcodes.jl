"""
`setup_barcodes(format::Symbol)`

Set up a Stateful iterator for the given barcode `format`. Accepts
- `haplotagging`
- `stlfr`
- `tellseq`
- `tenx`
"""
setup_barcodes(format::Symbol) = setup_barcodes(Val(format))

function setup_barcodes(::Val{:haplotagging})
    return Base.Iterators.Stateful(
        Base.Iterators.product('A', 1:96, 'C', 1:96, 'B', 1:96, 'D', 1:96)
    )
end

function setup_barcodes(::Val{:stlfr})
    return Base.Iterators.Stateful(
        Base.Iterators.product(1:2000, 1:2000, 1:2000)
    )
end

function setup_barcodes(::Val{:tellseq})
    return Base.Iterators.Stateful(
        Base.Iterators.product(("ATGC" for i in 1:18)...)
    )
end

function setup_barcodes(::Val{:tenx})
    return Base.Iterators.Stateful(
        Base.Iterators.product(("ATGC" for i in 1:16)...)
    )
end

setup_barcodes(::Val{:tenx}) = setup_barcodes(:tellseq)


"""
`get_next!(bc::T) where T<:Base.Iterators.Stateful`

Iterates to the next value of the barcode generator `bc` and formats
it into a proper string.
"""
function get_next!(bc::T)::String where T<:Base.Iterators.Stateful
    format_barcode(iterate(bc))
end

"""
`format_barcode(bc::Tuple{Tuple{Char, Int64, Char, Int64, Char, Int64, Char, Int64}, Nothing})`

Format a haplotagging style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
function format_barcode(bc::Tuple{Tuple{Char, Int64, Char, Int64, Char, Int64, Char, Int64}, Nothing})::String
    @inbounds bc[1][1] * lpad(bc[1][2],2, '0') * bc[1][3] * lpad(bc[1][4],2, '0') * bc[1][5] * lpad(bc[1][6],2, '0') * bc[1][7] * lpad(bc[1][8],2, '0')
    #=
    mapreduce(*, bc[1]) do x
        if x isa Int
            lpad(x, 2, "0")
        else
            x
        end
    end
    =#
end

"""
`format_barcode(bc::Tuple{Tuple{Int64, Int64, Int64}, Nothing})`

Format an stLFR style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
format_barcode(bc::Tuple{Tuple{Int64, Int64, Int64}, Nothing})::String = join(bc[1], "_")


"""
`format_barcode(bc::Tuple{NTuple{18, Char}, Nothing})`

Format a TELLseq/10X style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
format_barcode(bc::Tuple{NTuple{18, Char}, Nothing})::String = join(bc[1], "")
format_barcode(bc::Tuple{NTuple{16, Char}, Nothing})::String = join(bc[1], "")
format_barcode(bc::Nothing) = error("There are no more barcodes available for the selected barcode type.")
