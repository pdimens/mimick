"""
    setup_barcodes(format::Symbol) -> Iterators.Stateful{Iterators.product}

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


"""
    get_next!(bc::T) where T<:Base.Iterators.Stateful -> String

Iterates to the next value of the barcode generator `bc` and formats
it into a proper string.
"""
function get_next!(bc::T)::String where T<:Base.Iterators.Stateful
    format_barcode(iterate(bc))
end

"""
    format_barcode(bc::Tuple{Tuple{Char, Int64, Char, Int64, Char, Int64, Char, Int64}, Nothing}) -> String

Format a haplotagging style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
function format_barcode(bc::Tuple{Tuple{Char, Int64, Char, Int64, Char, Int64, Char, Int64}, Nothing})::String
    _bc = bc[begin]
    @inbounds _bc[1] * lpad(_bc[2],2, '0') * _bc[3] * lpad(_bc[4],2, '0') * _bc[5] * lpad(_bc[6],2, '0') * _bc[7] * lpad(_bc[8],2, '0')
end

"""
    format_barcode(bc::Tuple{Tuple{Int64, Int64, Int64}, Nothing}) -> String

Format an stLFR style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
format_barcode(bc::Tuple{Tuple{Int64, Int64, Int64}, Nothing})::String = join(bc[begin], "_")


"""
    format_barcode(bc::Tuple{NTuple{18, Char}, Nothing}) -> String

Format a TELLseq/10X style barcode for output. The input `bc` is expected to be created by iterating the barcode
generator from `setup_barcodes()`. Returns a `String`.
"""
format_barcode(bc::Tuple{NTuple{18, Char}, Nothing})::String = join(bc[begin], "")
format_barcode(bc::Tuple{NTuple{16, Char}, Nothing})::String = join(bc[begin], "")
format_barcode(bc::Nothing) = error("There are no more barcodes available for the selected barcode type.")
