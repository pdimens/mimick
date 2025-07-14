
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
get_sample_variants(vcf::String, sample_number::Int)

For a given `sample_number`, returns a Dict of the variants as a Vector of Pairs,
given as `position::Int => variant::String` e.g. `"chr1_1" => [192=>"T", 2427=>"A", 2582=>"C"]`
"""
function get_sample_variants(vcf::String, sample_number::Int)::Dict{String, Vector{Pair{Int64, String}}}
    contigs,ploidy = get_ploidy_and_contigs(vcf)
    out_dict = Dict{String, Vector{Pair{Int, String}}}()
    for contig in contigs
        for haplotype in 1:ploidy
            out_dict["$(contig)_$haplotype"] = Pair{Int, String}[]
        end
    end

    reader = VCF.Reader(open(vcf, "r"))
    @inbounds for record in reader
        @inbounds alt = VCF.alt(record)[1]
        chrom = VCF.chrom(record)
        pos = VCF.pos(record)
        gt = VCF.genotype(record, sample_number, "GT")
        alleles = parse.(UInt8, split(gt, ['/', '|']))
        for i in findall(!iszero, alleles)
            push!(out_dict["$(chrom)_$i"], (pos => alt))
        end
    end
    close(reader)
    return out_dict
end


#TODO there needs to be a way to dispatch on the mutation type, along with a universal way of encoding them

"""
mutate!()

Does an in-place mutation of `dna` given a mutation given as a `position::Int => nucleotide::Char`
"""
mutate!(dna::LongDNA{4}, mutation::Pair{Int,Char}) = @inbounds dna[mutation.first] = mutation.second

"""
mutate!(dna::LongDNA{4}, mutation::Pair{Int,String})

Does an in-place mutation of `dna` by deleting the nucleotides at `position::Int` to `position + length(nucleotides)`
"""
function mutate!(dna::LongDNA{4}, mutation::Pair{Int,String})
    deletionrange = mutation.first:mutation.first+length(mutation.second)
    deleteat!(dna, deletionrange)
    return dna
end
