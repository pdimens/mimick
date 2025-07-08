Random.seed!(123)

function find_nonoverlapping_ranges!(molecule::ProcessedMolecule, number_range::UnitRange{Int}, lengths::Vector{Int}, max_attempts::Int)
    range_start = number_range.start
    range_end = number_range.stop
    range_size = range_end - range_start + 1
    if sum(lengths) > range_size
        error("Cannot fit ranges of total length $(sum(lengths)) in range of size $range_size")
    end
    occupied_intervals = Vector{UnitRange{Int}}(undef, length(lengths))
    @inbounds for (i, len) in enumerate(lengths)
        placed = false
        attempts = 0
        last_possible = range_end - len

        while !placed && attempts < max_attempts
            attempts += 1
            start_pos = rand(range_start:last_possible)
            interval_attempt = start_pos:start_pos+len
            for _interval in occupied_intervals
                !isempty(intersect(interval_attempt,_interval)) && continue
            end
            occupied_intervals[i] = interval_attempt
            placed = true
        end
        if !placed
            println("Couldnt place $i")
            return Vector{UnitRange{Int}}(undef, 1)
            error("Could not place range of length $len after $max_attempts attempts")
        end
    end
    @inbounds molecule.read_breakpoints .= occupied_intervals
    return #occupied_intervals
end

"""
Given a `sequence`, a vector of `breakpoints` (as ranges, provided in the ProcessedMolecule), and read 1/2 lengths,
returns updates the Pair of Vectors in the ProcessedMolecule with the read sequences, accounting for circular.
"""
function extract_sequences!(molecule::ProcessedMolecule, sequence::LongDNA{4}, r1_len::Int, r2_len::Int)
    #reads = Vector{Pair{LongDNA{4}, LongDNA{4}}}(undef, length(molecule.read_breakpoints))
    #reads = Pair{Vector{LongDNA{4}}, Vector{LongDNA{4}}}(Vector{LongDNA{4}}(undef, molecule.read_breakpoints) => Vector{LongDNA{4}}(undef, molecule.read_breakpoints))
    r1_len -= 1
    r2_len -= 1
    @inbounds for (i,breakpoint) in enumerate(molecule.read_breakpoints)
        R1 = circular_index(sequence, breakpoint.start:breakpoint.start+r1_len)
        R2 = reverse(circular_index(sequence, (breakpoint.stop-r2_len):breakpoint.stop))
        if count(==(DNA_N), R1)/r1_len > 0.05 || count(==(DNA_N), R2)/r2_len > 0.05
            return Vector{Pair{LongDNA{4}, LongDNA{4}}}(undef, 1)
            error("Too many Ns")
        end
        @inbounds molecule.read_sequences.first[i] = R1
        @inbounds molecule.read_sequences.second[i] = R2
    end
    return # reads
end

"""
Circular indexing to account for end position overflow. Returns a LongDNA string
"""
@inline function circular_index(seq::LongDNA{4}, range::UnitRange{Int})::LongDNA{4}
    n = length(seq)
    if range.stop > length(seq)
        return seq[range.start:end] * seq[begin:(range.stop % n)]
    else
        return seq[range]
    end
end

"""
Given the length of a contig/interval and simulation parameters, randomly draws
a molecule length from an exponential distribution and returns it.
"""
function get_molecule_size(params::SimParams, seq_len::Int)::Int
    molecule_length = seq_len + 1
    # make sure the molecule is greater than 600 and less than the sequence
    while molecule_length > seq_len
        molecule_length = rand(params.molecule_length)
    end
    return trunc(Int, molecule_length)
end

"""
Given simulation parameters and molecule size, calculates how many how many fragments
and their lengths should be simulated from the molecule. Returns a Vector{Int} of 
fragment lengths.
"""
function calculate_fragments(params::SimParams, molsize::Int)::Vector{Int}
    # if a singleton proportion is provided, conditionally drop the number of reads to 1
    if params.singletons > 0.0 && rand() <= params.singletons
        N = 1.0
    else
        if params.molecule_coverage < 1
            # set a minimum number of reads to 2 to avoid singletons
            _n = max(2.0, ((molsize * params.molecule_coverage)) / (params.length_R1 + params.length_R2))
        else
            _n = params.molecule_coverage
        end
        # draw N from an exponential distribution with a minimum set to 2 reads to avoid singletons
        # set ceiling to avoid N being greater than can be sampled
        _exp = truncated(Exponential(_n), lower = 2.0, upper = molsize/(params.length_R1 + params.length_R2))
        N = max(2,rand(_exp))
    end
    n_reads = trunc(Int, N)
    return trunc.(Int,rand(params.insert_size, n_reads))
end

"""
Given a schema and parameters along with molecule length and number of reads, randomly generates
breakpoints for the molecule and the reads therein. Returns a `ProcessedMolecule`.
"""
function get_sequences(schema::Schema, params::SimParams, barcode::String, molecule_length::Int, fragments::Vector{Int})::ProcessedMolecule
    n_frags = length(fragments)
    molecule = ProcessedMolecule(schema, barcode, n_frags)
    molecule.barcode = barcode
    seq_len = length(schema.sequence)
    i = 0
    #seqs = Vector{Pair{LongDNA{4}, LongDNA{4}}}(undef, length(fragments))
    # set the max start position to be (length - mol_length) to avoid overflow if not circular
    adjusted_end = params.circular ? seq_len : seq_len - molecule_length

    while i <= params.attempts && !isassigned(molecule.read_sequences.first, n_frags)
        i += 1
        start_pos = rand(1:adjusted_end)
        end_pos = start_pos + molecule_length
        molecule.position = start_pos:end_pos
        # 50 attempts to create N reads with <95% ambiguous bases
        find_nonoverlapping_ranges!(molecule, molecule.position, fragments, 50)
        extract_sequences!(molecule, schema.sequence, params.length_R1, params.length_R2)
    end
    if !isassigned(molecule.read_sequences.first, n_frags)
        error("After $attempts attempts, unable to create reads for $(schema.chrom) from molecule spanning $start_pos-$end_pos. This could be due to either a failure to create $N non-overlapping reads or too many N's (ambiguous bases) in the resulting reads.")
    end
    #molecule.position = start_pos:end_pos,seqs
    return molecule
end
