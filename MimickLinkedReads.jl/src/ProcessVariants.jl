"""
    get_ploidy(vcf::String, :VCF) ->Tuple(contigs, ploidy)
    get_ploidy(bcf::String, :BCF) ->Tuple(contigs, ploidy)

Parse a BCF/VCf file and determine the ploidy from the genotype of the first sample at the first variant.
Returns a 2-Tuple of a Vector of contig names and the ploidy.
"""
get_ploidy_and_contigs(vcf::String, fmt::Symbol) = get_ploidy_and_contigs(vcf, Val(fmt))

function get_ploidy_and_contigs(vcf::String, ::Val{:VCF})::Tuple{Vector{String}, Int}
    reader = VCF.Reader(safe_read(vcf)) do reader
        _header = header(reader)
        X = findall(i -> metainfotag(_header.metainfo[i]) == "contig", 1:length(_header))
        contigs = [_header.metainfo[i]["ID"] for i in X]     
        for record in reader
            gt = VCF.genotype(record, 1, "GT")
            return (contigs, length(split(gt, ['/', '|'])))
        end
    end
end

function get_ploidy_and_contigs(bcf::String, ::Val{:BCF})::Tuple{Vector{String}, Int}
    reader = BCF.Reader(open(bcf, "r")) do reader
        _header = header(reader)
        X = findall(i -> metainfotag(_header.metainfo[i]) == "contig", 1:length(_header))
        contigs = [_header.metainfo[i]["ID"] for i in X]     
        for record in reader
            gt = BCF.genotype(record, 1, "GT")
            return (contigs, length(split(gt, ['/', '|'])))
        end
    end
end

"""
    get_samples(vcf::String, :VCF) -> Vector{String}
    get_samples(bcf::String, :BCF) -> Vector{String}

Parse a BCF/VCf file and determine the samples in it from the header
"""
get_samples(vcf::String, fmt::Symbol) = get_samples(vcf::String, Val(fmt))

function get_samples(vcf::String, ::Val{:VCF})::Vector{String}
    reader = VCF.Reader(safe_read(vcf)) do reader
        return header(reader).sampleID
    end
end

function get_samples(bcf::String, ::Val{:BCF})::Vector{String}
    reader = BCF.Reader(open(bcf, "r")) do reader
        return header(reader).sampleID
    end
end

"""
    vartype(ref::String, alt::Vector{String}) -> Symbol

Assess what type of variant (`:snp` or `:indel`) a locus
is based on what the `ref` and `alt` alleles are. Returns
a `Symbol`.
"""
function vartype(ref::String, alt::Vector{String})::Symbol
    l_ref = length(ref)
    l_alts = length.(alt)
    if l_ref == 1 && all(==(1), l_alts)
        return :snp
    else
        return :indel
    end
end

"""
    get_sample_variants(vcf::String, :VCF, sample_number::Int)
    get_sample_variants(bcf::String, :BCF, sample_number::Int)

For a given `sample_number`, returns a Dict of the variants as a Vector of Pairs,
given as `position::Int => variant::String` e.g. `"chr1_1" => [192=>"T", 2427=>"A", 2582=>"C"]`
"""
get_sample_variants(vcf::String, fmt::Symbol, sample_number::Int) = get_sample_variants(vcf, Val(fmt), sample_number)

function get_sample_variants(vcf::String, ::Val{:VCF}, sample_number::Int)::Dict{String, Vector{Mutation}}
    if !isfile(vcf)
        error("$vcf does not exist.")
    end
    contigs,ploidy = get_ploidy_and_contigs(vcf, :VCF)
    out_dict = Dict{String, Vector{Mutation}}()
    for contig in contigs
        @simd for haplotype in 1:ploidy
            out_dict["$(contig)_$haplotype"] = Mutation[]
        end
    end
    reader = VCF.Reader(safe_read(vcf))
    @inbounds for record in reader
        alt = VCF.alt(record)
        ref = VCF.ref(record)
        _vartype = vartype(ref, alt)
        chrom = VCF.chrom(record)
        pos = VCF.pos(record)
        gt = VCF.genotype(record, sample_number, "GT")
        alleles = parse.(UInt8, split(gt, ['/', '|']))
        @inbounds @simd for i in findall(!iszero, alleles)
            @inbounds push!(out_dict["$(chrom)_$i"], Mutation(pos, _vartype, ref, alt[alleles[i]]))
        end
    end
    close(reader)
    return out_dict
end

function get_sample_variants(bcf::String, ::Val{:BCF}, sample_number::Int)::Dict{String, Vector{Mutation}}
    if !isfile(bcf)
        error("$bcf does not exist.")
    end
    contigs,ploidy = get_ploidy_and_contigs(bcf, :BCF)
    out_dict = Dict{String, Vector{Mutation}}()
    for contig in contigs
        @simd for haplotype in 1:ploidy
            out_dict["$(contig)_$haplotype"] = Mutation[]
        end
    end
    reader = BCF.Reader(open(bcf, "r"))
    @inbounds for record in reader
        alt = BCF.alt(record)
        ref = BCF.ref(record)
        _vartype = vartype(ref, alt)
        chrom = BCF.chrom(record)
        pos = BCF.pos(record)
        gt = BCF.genotype(record, sample_number, "GT")
        alleles = parse.(UInt8, split(gt, ['/', '|']))
        @inbounds @simd for i in findall(!iszero, alleles)
            @inbounds push!(out_dict["$(chrom)_$i"], Mutation(pos, _vartype, ref, alt[alleles[i]]))
        end
    end
    close(reader)
    return out_dict
end

"""
    mutate!(seq, mutation::Mutation)
    mutate!(seq::LongDNA{4}, ::Val{:snp|:indel}, position::Int, replacement::String)

Does an in-place mutation of `seq` given a `Mutation`, which contains the position, type, and variant nucleotides.
Dispatches on different mutation types: `:snp`, `:indel`
"""
function mutate!(seq::LongDNA{4}, mutation::Mutation)
    @inbounds mutate!(seq, Val(mutation.type), mutation.position, mutation.ref, mutation.alt)
end

function mutate!(seq::LongDNA{4}, ::Val{:snp}, position::Int, ref::String, replacement::String)
    @inbounds seq[position] = replacement[1]
    return seq
end

function mutate!(seq::LongDNA{4}, ::Val{:indel}, position::Int, ref::String, replacement::String)
    position_end = position + length(ref) - 1 
    spliceinto!(seq, position:position_end, replacement)
end


"""
    spliceinto!(seq::BioSequence, span::UnitRange, x)
Delete the symbols at indices `span` in `seq`, and then copy `x` into the
first deleted position, then return `seq`.

`span` must be nonempty, or this function will throw an `ArgumentError`. To handle
potentially empty spans, check if the span is empty, and if so use `spliceinto(seq, first(span), x)`.

This was taken from BioSequences.jl, will be removed when the package is updated.
"""
function spliceinto!(seq::BioSequence, span::UnitRange, x)
    isempty(span) && throw(ArgumentError("span cannot be empty"))
    @boundscheck checkbounds(seq, span)
    oldlen = length(seq)
    xlen = length(x)
    if length(span) == xlen
        # Same lengths: Just copy in x
        copyto!(seq, first(span), x, 1, length(span))
    elseif length(span) < xlen
        # x is longer. Resize and shift to make room for more symbols,
        # then copy in x
        resize!(seq, oldlen + xlen - length(span))
        copyto!(seq, first(span) + xlen, seq, last(span) + 1, oldlen - last(span))
        copyto!(seq, first(span), x, 1, xlen)
    else
        # Span is longer. Delete the rightmost bases (to cause the smallest possible shift),
        # then copy in
        deleteat!(seq, first(span) + xlen:last(span))
        copyto!(seq, first(span), x, 1, xlen)
    end
    return seq
end

function build_sample_schema(schema::Dict{String, Schema}, variants::Dict{String, Vector{Mutation}})
    haplotypes = collect(keys(variants))
    contigs = collect(keys(schema))
    contigs_haplos = map(haplotypes) do i
        _x = findfirst(j -> occursin(j, i[begin:end-2]), contigs)
        (contigs[_x] => i)
    end
    d = Dict{String, Schema}()
    for contig in contigs_haplos
        _schema = schema[contig.first]
        _variants = variants[contig.second]
        _seq = copy(_schema.sequence)
        for variant in filter(i -> i.type == :snp, _variants)
            mutate!(_seq, variant)
        end
        # add the indels, except start from the end and work backwards so indels don't mess up downstream indices
        for variant in reverse(filter(i -> i.type != :snp, _variants))
            mutate!(_seq, variant)
        end
        d[contig.second] = Schema(
            parse(Int, split(contig.second, "_")[end]),
            _schema.chrom,
            _schema.tracker.reads_required,
            _seq
        )
    end
    return d
end