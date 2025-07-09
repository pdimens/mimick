struct BarcodeManifest
    barcodes
    output_type::String
    max::Int64
    function BarcodeManifest(barcodes, output_type::String, max::Int)
        new(Base.Iterators.Stateful(barcodes), setup_barcode_output(output_type), output_type, max)
    end
end

struct SimParams
    output_dir::String
    prefix::String
    error::Float32
    insert_size::Distribution
    insert_size_buffer::Int
    length_R1::Int
    length_R2::Int
    n_molecules::Union{UnitRange{Int},Distribution}
    molecule_length::Distribution
    molecule_coverage::Float64
    singletons::Float64
    attempts::Int
    circular::Bool
    function SimParams(prefix, error, read_distance, distance_stdev, length_R1, length_R2, n_molecules, molecule_length, molecule_coverage, singletons; circular::Bool = false, attempts::Int = 50)
        _prefix = Base.Filesystem.basename(prefix)
        _outdir = Base.Filesystem.dirname(prefix)
        if molecule_coverage < 1
            lowerbound = max(3000, trunc(Int, molecule_length * molecule_coverage * 0.25))
        else
            lowerbound = molecule_coverage * (read_distance + (0.75 * distance_stdev))
        end
        _exp = truncated(Exponential(molecule_length), lower = lowerbound)
        if n_molecules > 0
            _nmol = truncated(Exponential(n_molecules), lower = 1, upper = 2.5*n_molecules)
        else
            _nmol = abs(n_molecules):abs(n_molecules)
        end
        _ins = truncated(Normal(read_distance, distance_stdev), lower = 200)
        return new(_outdir, _prefix, error, _ins, read_distance, length_R1, length_R2, _nmol, _exp, molecule_coverage, singletons, attempts, circular)
    end
end

mutable struct SchemaTracker
    reads_current::Atomic{Int64}
    reads_required::Int
end

function Base.show(io::IO, data::SchemaTracker)
    println(io, "SchemaTracker Object")
    println(io, "  Reads Required: ", data.reads_required)
    println(io, "  Reads Current: ", data.reads_current)
end

"""
A struct storing all the necessary information to simulate reads from a contig or contig interval
"""
struct Schema
    haplotype::Int
    chrom::String
    tracker::SchemaTracker
    sequence::LongSequence{DNAAlphabet{4}}
    function Schema(haplotype::Int, chrom::String, reads_req::Int, sequence::LongSequence{DNAAlphabet{4}})
        new(haplotype, chrom, SchemaTracker(Atomic{Int64}(0),reads_req), sequence)
    end
end

function Base.show(io::IO, data::Schema)
    println(io, "Simulation Schema")
    for i in fieldnames(Schema)
        val = getfield(data, i)
        if i == :tracker
            println(io, "  tracker::SchemaTracker")
            println(io, "    reads_required::Int64 ", data.tracker.reads_required)
            println(io, "    reads_current::Int64 ", data.tracker.reads_current.value)
        elseif i == :sequence
            println(io, "  $i::", typeof(val), " length $(length(val))")
        else
            println(io, "  $i::", typeof(val), " ", val)
        end
    end
end

mutable struct ProcessedMolecule
    haplotype::Int
    barcode::String
    chrom::String
    position::UnitRange{Int}
    read_breakpoints::Vector{UnitRange{Int}}
    read_sequences::Pair{Vector{LongDNA{4}}, Vector{LongDNA{4}}}
end

function ProcessedMolecule(schema::Schema, barcode::String, n_reads::Int)
    ProcessedMolecule(schema.haplotype, barcode, schema.chrom, 0:0, Vector{UnitRange{Int}}(undef, n_reads), Pair{Vector{LongDNA{4}}, Vector{LongDNA{4}}}(Vector{LongDNA{4}}(undef, n_reads), Vector{LongDNA{4}}(undef, n_reads)))
end

function Base.show(io::IO, data::Union{ProcessedMolecule,SimParams})
    println(io, "$(typeof(data)) Object")
    for i in fieldnames(typeof(data))
        val = getfield(data, i)
        if typeof(val) <: Distribution
            println(io, " $i::$(supertype(typeof(val)))")
        else
            println(io, " $i::", typeof(val), " ", val)
        end
    end
end
