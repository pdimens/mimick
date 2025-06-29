mutable struct SchemaTracker
    reads_current::Int
    reads_required::Int
end

function Base.show(io::IO, data::SchemaTracker)
    println(io, "Object of type SchemaTracker")
    println(io, " Reads Required: ", data.reads_required)
    println(io, " Reads Current: ", data.reads_current)
end

"""
A struct storing all the necessary information to simulate reads from a contig or contig interval
"""
struct Schema
    haplotype::Int8
    chrom::String
    start_position::Int64
    end_position::Int64
    read_length::Int
    read_pairs_per_mol::Int
    tracker::SchemaTracker
    #mol_length::Int
    #mol_coverage::Float64
    #singletons::Float64
    #is_circular::Bool
    sequence::LongSequence{DNAAlphabet{4}}
    #function Schema(haplotype::Int, chrom::String, start_position::Int, end_position::Int, read_length::Int, read_pairs_per_mol::Int, reads_req::Int, mol_length::Int, mol_coverage::Float64, singletons::Float64, is_circular::Bool, sequence::LongSequence{DNAAlphabet{4}})
    #    new(Int8(haplotype), chrom, start_position, end_position, read_length, read_pairs_per_mol, SchemaTracker(0,reads_req), mol_length, mol_coverage, singletons, is_circular, sequence)
    function Schema(haplotype::Int, chrom::String, start_position::Int, end_position::Int, read_length::Int, read_pairs_per_mol::Int, reads_req::Int, sequence::LongSequence{DNAAlphabet{4}})
        new(Int8(haplotype), chrom, start_position, end_position, read_length, read_pairs_per_mol, SchemaTracker(0,reads_req), sequence)
    end
end

function Base.show(io::IO, data::Schema)
    println("Object of type Schema")
    for i in fieldnames(Schema)
        val = getfield(data, i)
        if i == :tracker
            println(io, " tracker::SchemaTracker")
            println(io, "  reads_required::Int64 ", data.tracker.reads_required)
            println(io, "  reads_current::Int64 ", data.tracker.reads_current)
        elseif i == :sequence
            println(io, " $i::", typeof(val), " length $(length(val))")
        else
            println(io, " $i::", typeof(val), " ", val)
        end
    end
end

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
function process_fasta(fasta::String, haplotype::Int64, coverage::Float64, mol_cov::Float64, mol_len::Int64, read_len::Vector{Int}, singletons::Float64, circular::Bool)::Vector{Schema}
    fai = index_fasta(fasta)
    mean_readlen = sum(read_len) ÷ 2
    mean_reads_per = max(1, mol_cov < 1 ? (mol_cov*mol_len)÷(sum(read_len)) : mol_cov)
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
            Schema(haplotype, chrom, 1, end_position, mean_readlen, trunc(Int,mean_reads_per), reads_required, mol_len, mol_cov, singletons, circular, seq)
        end
    end
end

"""
Iterate through input `fasta_files` to return a `Dict` of `Schema`
"""
function setup_schema(fasta_files::Vector{String}, coverage::Float64, mol_cov::Float64, mol_len::Int64, read_len::Vector{Int}, singletons::Float64, circular::Bool)::Dict{String, Schema}
    d = Dict{String, Schema}()
    for (i,j) in enumerate(fasta_files)
        schemas = process_fasta(j, i, coverage, mol_cov, mol_len, read_len, singletons, circular)
        for _schema  in schemas
            d[join([_schema.chrom,_schema.haplotype], "_")] = _schema
        end
    end
    return d
end
