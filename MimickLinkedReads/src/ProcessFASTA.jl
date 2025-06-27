mutable struct Schema
    haplotype::Int8
    chrom::String
    start_position::Int64
    end_position::Int64
    read_length::Int32
    read_pairs_per_mol::Int32
    reads_current::Int64
    reads_required::Int32
    mol_length::Int64
    mol_coverage::Float64
    singletons::Float32
    is_circular::Bool
    sequence::LongSequence{DNAAlphabet{4}}
end

function Base.show(io::IO, data::Schema)
    println("Object of type Schema")
    for i in fieldnames(Schema)
        val = getfield(data, i)
        println(io, " $i::", typeof(val), " ", val)
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

function process_fasta(fasta::String, haplotype::Int64, coverage::Float64, mol_cov::Float64, mol_len::Int64, read_len::Vector{Int}, singletons::Float64, circular::Bool)::Vector{Schema}
    fai = index_fasta(fasta)
    mean_readlen = sum(read_len) / 2
    mean_reads_per = max(1, mol_cov < 1 ? (mol_cov*mol_len)/(sum(read_len)) : mol_cov)
    FASTAReader(safe_open(fasta,"r"), index = fai, copy = false) do _fasta
        map(_fasta) do contig
            chrom = identifier(contig)
            seq = sequence(LongDNA{4}, contig)
            end_position = length(seq)
            normalized_length = end_position - count(==(DNA_N), seq)
            if normalized_length < 650
                println("Error in $fasta contig $chrom: contigs must have at least 650 unambiguous (non-N) bases.")
                #return
            end
            reads_required = trunc(Int,(coverage*normalized_length/mean_readlen)/2)
            Schema(haplotype, chrom, 1, end_position, mean_readlen, mean_reads_per, 0, reads_required, mol_len, mol_cov, singletons, circular, seq)
        end
    end
end

function setup_schema(fasta_files::Vector{String})::Dict
    d = Dict()
    for (i,j) in enumerate(fasta_files)
        schemas = process_fasta(j, i)
        for _schema  in schemas
            d[_schema.chrom * _schema.haplotype] = [_schema, 0]
        end
    end
    return d
end

# 

#def FASTAtoInventory(fasta, coverage, mol_cov, mol_len, read_len, singletons, circular) -> dict:
#    '''
#    Read the FASTA files and derive the contig name, start, and end positions and other simulation schema
#    and return a dict of Schema objects that's an inventory tracker in the form of
#    d[idx] = Schema
#    '''
#    inventory = {}
#    mean_reads_per = (mol_cov*mol_len)/(read_len*2) if mol_cov < 1 else mol_cov
#    mean_reads_per = max(1, mean_reads_per)
#    haplotype = 0
#    idx = 0
#    for _fasta in fasta:
#        haplotype += 1
#        with pysam.FastxFile(_fasta) as fa:
#            for contig in fa:
#                chrom = contig.name
#                start = 1
#                end = len(contig.sequence)
#                normalized_length = end - contig.sequence.count('N')
#                if normalized_length < 650:
#                    error_terminate(f"Error in {os.path.basename(fasta)} [yellow]contig {chrom}[/]: contigs must have at least 650 non-ambiguous ([yellow]N[/]) bases.")
#                reads_required = int((coverage*normalized_length/read_len)/2)
#                inventory[idx] = Schema(haplotype, chrom,start,end,read_len, mean_reads_per, reads_required, mol_len, mol_cov, singletons, circular, contig.sequence)
#                idx += 1
#    return inventory
