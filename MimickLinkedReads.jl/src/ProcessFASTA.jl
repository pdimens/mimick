"""
Delete an existing faidx (if present) and re-index the input fasta file.
Returns the name of the fasta index file.
"""
function index_fasta(fasta::String)::String
    fai = fasta * ".fai"
    if isfile(fai)
        rm(fai)
    end
    #TODO DEAL WITH GZIP SITUATION
    try
        faidx(fasta)
    catch y
        println("Failed to index $fasta, is it a properly formatted FASTA file?")
    finally
        return fai
    end
end

"""
Read an input fasta and return a Schema object.
"""
function process_fasta(fasta::String, haplotype::Int64, coverage::Float64, read_len::Vector{Int})::Vector{Schema}
    fai = index_fasta(fasta)
    mean_readlen = sum(read_len) ÷ 2
    FASTAReader(safe_read(fasta), index = fai, copy = false) do _fasta
        map(_fasta) do contig
            chrom = String(identifier(contig))
            seq = sequence(LongDNA{4}, contig)
            end_position = length(seq)
            normalized_length = end_position - count(==(DNA_N), seq)
            if normalized_length < 600
                println("Error in $fasta contig $chrom: contigs must have at least 600 unambiguous (non-N) bases.")
                #return
            end
            reads_required = trunc(Int,(coverage*normalized_length/mean_readlen)/2)
            Schema(haplotype, chrom, reads_required, seq)
        end
    end
end

"""
Iterate through input `fasta_files` to return a `Dict` of `Schema`
"""
function setup_schema(fasta_files::Vector{String}, coverage::Float64, read_len::Vector{Int})::Dict{String, Schema}
    d = Dict{String, Schema}()
    @inbounds for (i,j) in enumerate(fasta_files)
        schemas = process_fasta(j, i, coverage, read_len)
        for _schema  in schemas
            d[join([_schema.chrom,_schema.haplotype], "_")] = _schema
        end
    end
    return d
end
