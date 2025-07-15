"""
`get_ploidy(vcf::String)`

Parse a VCf file and determine the ploidy from the genotype of the first sample at the first variant.
Returns the list of contig names and the ploidy.
"""
function get_ploidy_and_contigs(vcf::String)::Tuple{Vector{String}, Int}
    reader = VCF.Reader(open(vcf, "r")) do reader
        _header = header(reader)
        X = findall(i -> metainfotag(_header.metainfo[i]) == "contig", 1:length(_header))
        contigs = [_header.metainfo[i]["ID"] for i in X]     
        for record in reader
            gt = VCF.genotype(record, 1, "GT")
            return (contigs, length(split(gt, ['/', '|'])))
        end
    end
end

"""
`get_samples(vcf::String)`

Parse a VCf file and determine the samples in it from the header
"""
function get_samples(vcf::String)::Vector{String}
    reader = VCF.Reader(open(vcf, "r")) do reader
        return header(reader).sampleID
    end
end

"""
vartype(ref::String, alt::Vector{String})

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
get_sample_variants(vcf::String, sample_number::Int)

For a given `sample_number`, returns a Dict of the variants as a Vector of Pairs,
given as `position::Int => variant::String` e.g. `"chr1_1" => [192=>"T", 2427=>"A", 2582=>"C"]`
"""
function get_sample_variants(vcf::String, sample_number::Int)::Dict{String, Vector{Mutation}}
    contigs,ploidy = get_ploidy_and_contigs(vcf)
    out_dict = Dict{String, Vector{Mutation}}()
    for contig in contigs
        @simd for haplotype in 1:ploidy
            out_dict["$(contig)_$haplotype"] = Mutation[]
        end
    end
    reader = VCF.Reader(open(vcf, "r"))
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

#TODO THIS DOESNT WORK FOR SOME REASON
function mutate!(seq::LongDNA{4}, ::Val{:indel}, position::Int, ref::String, replacement::String)
    # remove the ref base(s)
    deleteat!(seq, position:length(ref)-1)
    # add the alt base(s) one by one bc a built-in method doesn't yet exist
    for (pos,i) in enumerate(replacement)
        insert!(seq, position + (pos-1), i)
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
    #return contigs_haplos
    d = Dict{String, Schema}()
    for contig in contigs_haplos
        _schema = schema[contig.first]
        _variants = variants[contig.second]
        _seq = copy(_schema.sequence)
        for variant in filter(i -> i.type == :snp, _variants)
            mutate!(_seq, variant)
        end
        # add the indels, except start backwards to seek the indicies
        # from clashing with each other
        #= disabled for now b/c indels aren't mutating properly
        for variant in reverse(filter(i -> i.type != :snp, _variants))
            mutate!(_seq, variant)
        end
        =#
        d[contig.second] = Schema(
            parse(Int, split(contig.second, "_")[end]),
            _schema.chrom,
            _schema.tracker.reads_required,
            _seq
        )
    end
    return d
end