"""
`index_fasta(fasta::String)`

## DEPRECATED
Delete an existing faidx (if present) and re-index the input fasta file.
Returns the name of the fasta index file.
"""
function index_fasta(fasta::String)::String
    if !isfile(fasta)
        error("$fasta was not found.")
    end
    fai = fasta * ".fai"
    if isfile(fai)
        rm(fai)
    end
    try
        faidx(fasta)
    catch y
        error("Failed to index $fasta, is it a properly formatted FASTA file? B/Gzipped FASTA are not supported.")
    end
    return fai
end

"""
    mask_N!(seq::LongDNA{4})

Replace instances of N with a random base. Mutates `seq` in place.
"""
function mask_N!(seq::LongDNA{4})
    for i in eachindex(seq)
        if seq[i] == DNA_N
            seq[i] = rand(ACGT)
        end
    end
    return nothing
end


"""
`process_fasta(fasta::String, haplotype::Int, coverage::Float64, read_len::Vector{Int}; mask::Bool=False)`

Read an input fasta and return a Vector of `Schema` objects, one for each contig in the file.
Use `mask=true` to replace Ns with random bases.
"""
function process_fasta(fasta::String, haplotype::Int, coverage::Float64, read_len::Vector{Int}; mask::Bool = false)::Vector{Schema}
    #fai = index_fasta(fasta)
    mean_readlen = sum(read_len) ÷ 2
    FASTAReader(safe_read(fasta), index=nothing, copy=false) do _fasta
        map(_fasta) do contig
            chrom = String(identifier(contig))
            seq = sequence(LongDNA{4}, contig)
            if mask
                mask_N!(seq)
            end
            end_position = lastindex(seq)
            normalized_length = end_position - count(==(DNA_N), seq)
            if normalized_length < 600
                error("Error in $fasta contig $chrom: contigs must have at least 600 unambiguous (non-N) bases.")
            end
            reads_required = trunc(Int, (coverage * normalized_length / mean_readlen) / 2)
            Schema(haplotype, chrom, reads_required, seq)
        end
    end
end

"""
`setup_schema(fasta_files::Vector{String}, coverage::Union{Int,Float64}, read_len::Vector{Int}, mol_length::Int; mask::Bool=false)`
`setup_schema(fasta_file::String, coverage::Union{Int,Float64}, read_len::Vector{Int}, mol_length::Int; mask::Bool=false)`

Iterate through input `fasta_files` or `fasta_file` to return a `Dict` of `Schema`. If given
a single FASTA input, all Schema will have `.haplotype` set to `0`. If given a Vector of FASTA
files, will encode the haplotypes based on order provided, along with having the haplotype suffixing
the key with an underscore (e.g. `chr1_1 => Schema`, `chr1_2 => Schema`, etc.). `mol_length` is used
to check sequence length and output a notice about handling contigs smaller than specific mean molecule length.
Use `mask=true` to replace Ns in the fasta file(s) with random ATCG bases.
"""
setup_schema(fasta_file::String, coverage::Union{Int,Float64}, read_len::Vector{Int}, mol_length::Int; mask::Bool=false)::Dict{String,Schema} = setup_schema(String[fasta_file], Float64(coverage), read_len, mol_length, mask = mask)
setup_schema(fasta_files::Vector{String}, coverage::Int, read_len::Vector{Int}, mol_length::Int; mask::Bool=false)::Dict{String,Schema} = setup_schema(fasta_files, Float64(coverage), read_len, mol_length, mask=mask)
function setup_schema(fasta_files::Vector{String}, coverage::Union{Int,Float64}, read_len::Vector{Int}, mol_length::Int; mask::Bool=false)::Dict{String,Schema}
    d = Dict{String,Schema}()
    @inbounds for (i, j) in enumerate(fasta_files)
        schemas = process_fasta(j, i, Float64(coverage), read_len, mask = mask)
        @simd for _schema in schemas
            d[join([_schema.chrom, _schema.haplotype], "_")] = _schema
        end
    end
    for i in values(d)
        if i.sequence.len <= mol_length
            @info "Contigs with lengths <= $mol_length (mean molecule length) will likely result in molecules that span the entire contig."
            break
        end
    end
    return d
end
