"""
    mimick(fasta::Vector{String}, format::String; kwargs...) -> Nothing

The wrapper function that simulates linked-reads in a given `format` for a single sample using the input haplotypes provided as a Vector to `fasta`.
Writes one pair of R1 and R2 reads. The `format` is expected to be one of `"haplotagging"`, `"stlfr"`, `"tellseq"`, `"tenx"`, `"standard"`,
`"standard:haplotagging"`, `"standard:stlfr"`

## Keyword Arguments
- prefix: `String` of the output file prefix, can be in the form of `path/to/dir/pref` (default: `simulated/SIM`)
- coverage: `Int` or `Float` of desired output read depth (default: `30`)
- n_molecules: `Int` of the average number of molecules per barcode. A negative number will be fixed for that amount (e.g. `-1` = "exactly one molecules per barcode", default: `2`)
- mol_cov: `Float64` of proportional coverage (e.g. `0.3` = 30%) or an `Int` of average reads per molecule (default: `0.2`)
- mol_len: `Int` of average molecule length, in bp (default: `80000`)
- insert_length: `Int` of average length of inserts (molecule fragments) from which sequences are generated (default: `500`)
- insert_stdev: `Int` of standard deviation of `insert_length` (default: `50`)
- read_length: `Vector{Int}` of [R1,R2] read lengths in bp (default: `[150,150]`)
- singletons: `Float64` if the proportion of barcodes that will have only one read pair (default: `0.35`)
- circular: `Bool` of whether to treat input FASTA contigs as circular DNA when simulating molecules (default: `false`)
- attempts: `Int` of how many times to attempt creating a feasible molecule for a barcode before terminating with an error (default: `25`)
- seed: `Int` of seed for randomization. Numbers >=0 will set the seed to that value, whereas negative numbers ignore setting a seed (default: `-1`)
"""
function mimick(fasta::Vector{String}, format::String; prefix::String = "simulated/SIM", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_length::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = -1)
    if seed >= 0
        Random.seed!(seed)
    end
    if length(read_length) > 2
        @info "More than 2 read lengths were provided. Only using the first two for R1 and R2 reads, respectively"
    elseif length(read_length) == 1
        @info "Only 1 read lengths was provided. Using it for both R1 and R2 read lengths"
        read_length = [read_length, read_length]
    end
    schema = setup_schema(fasta, coverage, read_length)
    bc_fmt, fq_fmt = interperet_format(format)
    barcodes = setup_barcodes(bc_fmt)
    params = SimParams(prefix,insert_length, insert_stdev, read_length[1], read_length[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
    mkpath(dirname(prefix))
    n_haps = length(schema)
    pbar = ProgressBar(; transient = true)
    job = addjob!(pbar; N = n_haps, description="Simulating Reads")
    with(pbar) do
        open(GzipCompressorStream,"$prefix.R1.fq.gz", "w") do R1; open(GzipCompressorStream, "$prefix.R2.fq.gz", "w") do R2 
            while !isempty(schema)
                candidates = keys(schema)
                n_mol = get_n_molecules(params)
                molecule_targets = rand(candidates, n_mol)
                for target in molecule_targets
                    molsize = get_molecule_size(params, length(schema[target].sequence))
                    frags = calculate_insert_sizes(params, molsize)
                    molecule = get_sequences(schema[target], params, get_next!(barcodes), molsize, frags)
                    write(R1, format_R1(fq_fmt, molecule))
                    write(R2, format_R2(fq_fmt, molecule))
                    schema[target].tracker.reads_current += length(frags)
                end
                filter!(is_incomplete, schema)
                Progress.update!(job, i = n_haps - length(schema))
            end
        end;end
    end
end

"""
    mimick(fasta::String, vcf:String; kwargs...) -> Nothing

The wrapper function that simulates linked-reads in a given `format` for all the samples in `vcf`. There is a single `fasta` input file that each sample
will have converted into its haplotypes by applying the SNP and INDEL variants in `vcf`. Assumes phased genotypes. Writes one pair of R1 and R2 reads per sample.

The `format` is expected to be one of `"haplotagging"`, `"stlfr"`, `"tellseq"`, `"tenx"`, `"standard"`, `"standard:haplotagging"`, `"standard:stlfr"`

## Keyword Arguments
- prefix: `String` of the output file prefix, can be in the form of `path/to/dir/pref` (default: `simulated/SIM`)
- coverage: `Int` or `Float` of desired output read depth (default: `30`)
- n_molecules: `Int` of the average number of molecules per barcode. A negative number will be fixed for that amount (e.g. `-1` = "exactly one molecules per barcode", default: `2`)
- mol_cov: `Float64` of proportional coverage (e.g. `0.3` = 30%) or an `Int` of average reads per molecule (default: `0.2`)
- mol_len: `Int` of average molecule length, in bp (default: `80000`)
- insert_length: `Int` of average length of inserts (molecule fragments) from which sequences are generated (default: `500`)
- insert_stdev: `Int` of standard deviation of `insert_length` (default: `50`)
- read_length: `Vector{Int}` of [R1,R2] read lengths in bp (default: `[150,150]`)
- singletons: `Float64` if the proportion of barcodes that will have only one read pair (default: `0.35`)
- circular: `Bool` of whether to treat input FASTA contigs as circular DNA when simulating molecules (default: `false`)
- attempts: `Int` of how many times to attempt creating a feasible molecule for a barcode before terminating with an error (default: `25`)
- seed: `Int` of seed for randomization. Numbers >=0 will set the seed to that value, whereas negative numbers ignore setting a seed (default: `-1`)
"""
function mimick(fasta::String, vcf::String, format::String; prefix::String = "simulated/", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_length::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = -1)
    if seed >= 0
        Random.seed!(seed)
    end
    master_schema = setup_schema(fasta, coverage, read_length)
    params = SimParams(prefix, insert_length, insert_stdev, read_length[1], read_length[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
    bc_fmt, fq_fmt = interperet_format(format)
    mkpath(dirname(prefix))
    samplenames = get_samples(vcf)
    pbar = ProgressBar(; transient = true)
    job = addjob!(pbar; N=length(samplenames), description="Simulating Samples")
    with(pbar) do
        @inbounds @sync for idx in eachindex(samplenames)
            Base.Threads.@spawn begin
                outprefix = prefix * samplenames[idx]
                barcodes = setup_barcodes(bc_fmt)
                variants = get_sample_variants(vcf, idx)
                schema = build_sample_schema(master_schema, variants)
                # variants are no longer needed, free it up
                variants = nothing
                open(GzipCompressorStream, "$outprefix.R1.fq.gz", "w") do R1; open(GzipCompressorStream, "$outprefix.R2.fq.gz", "w") do R2 
                    while !isempty(schema)
                        candidates = keys(schema)
                        n_mol = get_n_molecules(params)
                        molecule_targets = rand(candidates, n_mol)
                        for target in molecule_targets
                            molsize = get_molecule_size(params, length(schema[target].sequence))
                            frags = calculate_insert_sizes(params, molsize)
                            molecule = get_sequences(schema[target], params, get_next!(barcodes), molsize, frags)
                            write(R1, format_R1(fq_fmt, molecule))
                            write(R2, format_R2(fq_fmt, molecule))
                            schema[target].tracker.reads_current += length(frags)
                        end
                        filter!(is_incomplete, schema)
                    end
                end;end
                Progress.update!(job)
            end
        end
    end
end
