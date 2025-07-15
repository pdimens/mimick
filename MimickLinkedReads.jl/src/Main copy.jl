using ArgParse
using ArgParse
using Base.Threads
using BioSequences
using CodecZlib
using Distributions
using FASTX
using Random

include("Structs.jl")
include("Barcodes.jl")
include("Common.jl")
include("FormatFASTQ.jl")
include("Breakpoints.jl")
include("ProcessFASTA.jl")
include("ProcessFASTQ.jl")

function parse_commandline()
    s = ArgParseSettings(
        "Simulate linked-read data",
        version = "2.0-beta",
    )
    add_arg_group(s, "Positional Arguments")
    @add_arg_table s begin
        "fasta"
        arg_type = String
        nargs = '+'
        required = true
        help = "haplotypes to simulate reads from, in FASTA format"
        range_tester = isfile
    end
    add_arg_group(s, "General Options")
    @add_arg_table s begin
        "--format", "-f"
            arg_type = String
            default = "haplotagging"
            help = "Linked-read output style (10x, haplotagging, stlfr, tellseq, standard, standard_haplotagging, standard_stlfr)"
        "--output-prefix", "-o"
            arg_type = String
            default = "simulated/SIM"
            help = "filename prefix for output files"
    end
    add_arg_group(s, "Read Simulation Parameters")
    @add_arg_table s begin
        "--coverage", "-c"
            arg_type = Float64
            default = 30
            help = "total sequencing depth for each haplotype"
        "--insert-length", "-i"
            arg_type = Int
            default = 500
            range_tester = >(0)
            help = "average insert length, in bp"
        "--insert-stdev", "-d"
            arg_type = Int
            default = 50
            help = "standard deviation of insert length"
        "--read-lengths", "-l"
            arg_type = Int
            nargs = 2
            default = [150,150]
            help = "the lengths of the R1 (forward) and R2 (reverse) reads, expects 2 numbers"
        "--seed", "-s"
            arg_type = Int
            default = 0
            help = "set a seed for the random number generators"
    end
    add_arg_group(s, "Linked Read Parameters")
    @add_arg_table s begin
        "--circular", "-O"
            help = "treat contigs in FASTA files as circular"
            action = :store_true
        "--molecules-per", "-N"
            arg_type = Int
            default = 2
            help = "if positive: average number of molecules per barcode (fixed for a number if negative)"
        "--molecule-coverage", "-C"
            arg_type = Float64
            default = 0.2
            help = "<1 is the proportion of molecule covered, >1 is an average of that many reads per molecule"
        "--molecule-length", "-L"
            arg_type=Int64
            default = 80000
            help = "average molecule length, in bp"
        "--singletons", "-S"
            arg_type = Float64
            default = 0.35
            help = "proportion of barcodes that will be singletons (unlinked)"
        "--molecule-attempts", "-A"
            arg_type = Int
            default = 25
            help = "number of attempts to randomly simulate a molecule before terminating"
    end

    return parse_args(s)
end

"""
The main function that takes arguments and does the thing
"""
function mimick(fasta_files::Vector{String}; format::String, prefix::String = "simulated/SIM", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = 0)
    if seed > 0
        Random.seed!(seed)
    end
    schema = setup_schema(fasta_files, coverage, read_len)
    barcodes = setup_barcodes(Symbol(format))
    params = SimParams(prefix, 0.001, insert_length, insert_stdev, read_len[1], read_len[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
    mkpath(dirname(prefix))
    writer = FastqWriter(prefix, Symbol(format))
    total_reads = 0
    for i in keys(schema)
        total_reads += schema[i].tracker.reads_total
    end
    pbar = ProgressBar()
    job = addjob!(pbar; N=total_reads)
    start!(pbar)
    
    while !isempty(schema)
        candidates = keys(schema)
        n_mol = get_n_molecules(params)
        molecule_targets = rand(candidates, n_mol)
        for target in molecule_targets
            molsize = get_molecule_size(params, length(schema[target].sequence))
            frags = calculate_fragments(params, molsize)
            molecule = get_sequences(schema[target], params, get_next!(barcodes), molsize, frags)
            Base.Threads.atomic_add!(schema[target].tracker.reads_current, length(frags))
            submit!(writer, molecule)
            update!(job, length(frags))
            render(pbar)
        end
        update!(schema)
    end
    stop!(writer)
    stop!(pbar)
end

function main()
    parsed_args = parse_commandline()
    mimick(
        parsed_args["fasta"],
        format = parsed_args["format"],
        prefix = parsed_args["output-prefix"],
        coverage = parsed_args["coverage"],
        n_molecules = parsed_args["molecules-per"],
        mol_cov = parsed_args["molecule-coverage"],
        mol_len = parsed_args["molecule-length"],
        insert_length = parsed_args["insert-length"],
        insert_stdev = parsed_args["insert-stdev"],
        read_len = parsed_args["read-lengths"],
        singletons = parsed_args["singletons"],
        circular = parsed_args["circular"],
        attempts = parsed_args["molecule-attempts"],
        seed = parsed_args["seed"]
    )    
end

main()

