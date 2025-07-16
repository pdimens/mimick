"""
    mimick(fasta_files::Vector{String}; format::String, prefix::String = "simulated/SIM", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = -1

The main function that takes arguments and does the thing. Setting a `seed` to `<0` (the default) is a shortcut not to use any seed at all.
"""
function mimick(fasta_files::Vector{String}; format::String, prefix::String = "simulated/SIM", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = 0)
    if seed >= 0
        Random.seed!(seed)
    end
    schema = setup_schema(fasta_files, coverage, read_len)
    bc_fmt, fq_fmt = interperet_format(format)
    barcodes = setup_barcodes(bc_fmt)
    params = SimParams(prefix,insert_length, insert_stdev, read_len[1], read_len[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
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
                    frags = calculate_fragments(params, molsize)
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
    mimick(fasta_file::String, vcf_file::String; format::String, prefix::String = "simulated/", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = 0)


The main function that takes arguments and does the thing, but for one fasta file and a VCF of variants.
"""
function mimick(fasta_file::String, vcf_file::String; format::String, prefix::String = "simulated/", coverage::Union{Int,Float64} = 30, n_molecules::Int = 2, mol_cov::Float64 = 0.2, mol_len::Int64 = 80000, insert_length::Int = 500, insert_stdev::Int = 50, read_len::Vector{Int} = [150,150], singletons::Float64 = 0.35, circular::Bool = false, attempts::Int = 25, seed::Int = -1)
    if seed >= 0
        Random.seed!(seed)
    end
    master_schema = setup_schema(fasta_file, coverage, read_len)
    params = SimParams(prefix, insert_length, insert_stdev, read_len[1], read_len[2], n_molecules, mol_len, mol_cov, singletons; circular = circular, attempts = attempts)
    bc_fmt, fq_fmt = interperet_format(format)
    mkpath(dirname(prefix))
    samplenames = get_samples(vcf_file)
    pbar = ProgressBar(; transient = true)
    job = addjob!(pbar; N=length(samplenames), description="Simulating Samples")
    with(pbar) do
        @inbounds @sync for idx in eachindex(samplenames)
            Base.Threads.@spawn begin
                outprefix = prefix * samplenames[idx]
                barcodes = setup_barcodes(bc_fmt)
                variants = get_sample_variants(vcf_file, idx)
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
                            frags = calculate_fragments(params, molsize)
                            molecule = get_sequences(schema[target], params, get_next!(barcodes), molsize, frags)
                            #for record in format_R1(format, molecule)
                            write(R1, format_R1(fq_fmt, molecule))
                            #end
                            # convert to fastq records and write R2
                            #for record in format_R2(format, molecule)
                            write(R2, format_R2(fq_fmt, molecule))
                            #end
                            #return molecule
                            schema[target].tracker.reads_current += length(frags)
                        end
                        filter!(is_incomplete, schema)
                        #update!(schema)
                    end
                end;end
                Progress.update!(job)
            end
        end
    end
end
