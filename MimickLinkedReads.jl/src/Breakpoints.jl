"""
`cumsum_to_ranges(values::Vector{Int64})`

Shortcut method of finding fragment breakpoints. The first value of the `values` Vector
is expected to be the molecule start point. The subsequent values (indicies 2+) are expected
to be the lengths of the fragments that need breakpoints. This method will return sequential
non-overlapping breakpoints as a `Vector{UnitRange{Int}}`

"""
function cumsum_to_ranges(values::Vector{Int})::Vector{UnitRange{Int}}
    cs = cumsum(values)
    @inbounds begin
        return [cs[i]+1:cs[i+1] for i in eachindex(cs)[begin:end-1]]
    end
end

"""
find_nonoverlapping_ranges!(molecule::ProcessedMolecule, lengths::Vector{Int}, max_attempts::Int)

Given the molecule breakpoints in the `ProcessedMolecule` and the `lengths` of the inserts, will
attempt (up to `max_attempts`) to randomly identify non-overlapping breakpoints. Mutates `molecule.breakpoints`
in place.
"""
function find_nonoverlapping_ranges!(molecule::ProcessedMolecule, lengths::Vector{Int}, max_attempts::Int)
    range_start = molecule.position.start
    range_end = molecule.position.stop
    range_size = range_end - range_start + 1
    
    # if the lengths are kinda close to the total molecule size, just use them sequentially and skip everything else
    if sum(lengths) > range_size - 200
        @inbounds molecule.read_breakpoints .= cumsum_to_ranges([range_start, lengths...])
        return
    end
    occupied_intervals = Vector{UnitRange{Int}}(undef, length(lengths))
    @inbounds for (i, len) in enumerate(lengths)
        placed = false
        attempts = 0
        last_possible = range_end - len

        while !placed && attempts < max_attempts
            attempts += 1
            start_pos = rand(range_start:last_possible)
            interval_attempt = start_pos:start_pos + len
            for _interval in occupied_intervals
                !isempty(intersect(interval_attempt,_interval)) && continue
            end
            @inbounds occupied_intervals[i] = interval_attempt
            placed = true
        end
        if !placed
            println("Couldnt place $i")
            return Vector{UnitRange{Int}}(undef, 1)
            error("Could not place range of length $len after $max_attempts attempts")
        end
    end
    @inbounds molecule.read_breakpoints .= occupied_intervals
    return
end

"""
extract_sequences!(molecule::ProcessedMolecule, sequence::LongDNA{4}, r1_len::Int, r2_len::Int)

Given a `sequence`, a `ProcessedMolecule` containing breakpoints and read 1/read 2 lengths,
mutates `ProcessedMolecule.read_sequences` in place with updated read sequences.
"""
function extract_sequences!(molecule::ProcessedMolecule, sequence::LongDNA{4}, r1_len::Int, r2_len::Int)
    r1_len -= 1
    r2_len -= 1
     @fastmath @inbounds for (i,breakpoint) in enumerate(molecule.read_breakpoints)
        R1 = circular_index(sequence, breakpoint.start:breakpoint.start+r1_len)
        R2 = reverse(circular_index(sequence, (breakpoint.stop-r2_len):breakpoint.stop))
        if count(==(DNA_N), R1)/r1_len > 0.05 || count(==(DNA_N), R2)/r2_len > 0.05
            return Vector{Pair{LongDNA{4}, LongDNA{4}}}(undef, 1)
            error("Too many Ns")
        end
        @inbounds molecule.read_sequences.first[i] = R1
        @inbounds molecule.read_sequences.second[i] = R2
    end
    return
end

"""
`circular_index(seq::LongDNA{4}, range::UnitRange{Int})`

Circular indexing to account for end position overflow. Returns a `LongDNA{4}` string
"""
@inline function circular_index(seq::LongDNA{4}, range::UnitRange{Int})::LongDNA{4}
    n = length(seq)
    if range.stop > n
        @inbounds return seq[range.start:end] * seq[begin:(range.stop % n)]
    else
        @inbounds return seq[range]
    end
end

"""
`get_molecule_size(params::SimParams, seq_len::Int)`

Given the length of a contig/interval and simulation parameters, randomly draws
a molecule length from an exponential distribution and returns it as an `Int`.
"""
function get_molecule_size(params::SimParams, seq_len::Int)::Int
    molecule_length = seq_len + 1
    # make sure the molecule is not longer than the sequence
    while molecule_length > seq_len
        molecule_length = rand(params.molecule_length)
    end
    return trunc(Int, molecule_length)
end

"""
`calculate_fragments(params::SimParams, molsize::Int)`

Given simulation parameters and molecule size, calculates how many how many fragments
and their lengths should be simulated from the molecule. Returns a `Vector{Int}` of 
fragment lengths.
"""
function calculate_fragments(params::SimParams, molsize::Int)::Vector{Int}
    # if a singleton proportion is provided, conditionally drop the number of reads to 1
    if params.singletons > 0.0 && rand() <= params.singletons
        N = 1.0
    else
        if params.molecule_coverage < 1
            # set a minimum number of 2 reads to avoid singletons
            @fastmath exact_inserts = (molsize * params.molecule_coverage) / params.insert_size_buffer
            _n = max(2.0, exact_inserts)
        else
            # molecule coverage given as an integer, convert to float
            _n = Float64(params.molecule_coverage)
        end
        # the number of fragments should be imperfect, so we draw N from an
        # exponential distribution with a minimum set to 2 reads to avoid singletons
        # set ceiling to avoid N being greater than can be sampled
        @fastmath max_possible = molsize/params.insert_size_buffer
        _exp = truncated(Exponential(_n), lower = 2.0, upper = max_possible)
        N = rand(_exp)
    end
    n_reads = trunc(Int, N)
    read_sizes = rand(params.insert_size, n_reads)
    attempts = 0
    while sum(read_sizes) > molsize && attempts < params.attempts
        attempts += 1
        read_sizes = rand(params.insert_size, n_reads)
    end
    if sum(read_sizes) > molsize
        error("FAILED TO FIND PRACTICAL READ SIZES LESS THAN MOLECULE.\nMolecule Length: $molsize\nRead Lengths: $read_sizes\nTotal Insert Length:$(sum(read_sizes))")
    end
    return trunc.(Int, read_sizes)
end

"""
`get_sequences(schema::Schema, params::SimParams, barcode::String, molecule_length::Int, fragments::Vector{Int})`

Given a schema and parameters along with molecule length and number of reads, randomly generates
breakpoints for the molecule and the reads therein. Returns a `ProcessedMolecule`, which contains
all the necessary information to convert the sequences into barcoded FASTQ format.
"""
function get_sequences(schema::Schema, params::SimParams, barcode::String, molecule_length::Int, fragments::Vector{Int})::ProcessedMolecule
    n_frags = length(fragments)
    molecule = ProcessedMolecule(schema, barcode, n_frags)
    molecule.barcode = barcode
    seq_len = length(schema.sequence)
    i = 0
    # set the max start position to be (length - mol_length) to avoid overflow if not circular
    adjusted_end = params.circular ? seq_len : seq_len - molecule_length

    while i <= params.attempts && !isassigned(molecule.read_sequences.first, n_frags)
        i += 1
        start_pos = rand(1:adjusted_end)
        end_pos = start_pos + molecule_length
        molecule.position = start_pos:end_pos
        # 50 attempts to create N reads with <95% ambiguous bases
        find_nonoverlapping_ranges!(molecule, fragments, 50)
        extract_sequences!(molecule, schema.sequence, params.length_R1, params.length_R2)
    end
    if !isassigned(molecule.read_sequences.first, n_frags)
        error("After $attempts attempts, unable to create reads for $(schema.chrom) from molecule spanning $start_pos-$end_pos. This could be due to either a failure to create $N non-overlapping reads or too many N's (ambiguous bases) in the resulting reads.")
    end
    return molecule
end
