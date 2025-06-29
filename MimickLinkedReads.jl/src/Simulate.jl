struct Molecule
    haplotype::Int8
    chrom::String
    start_position::Int
    end_position::Int
    length::Int
    barcode::String
    output_barcode::String
    read_count::Int
end

function find_breakpoints(schema::Schema, params::SimParams )
    molecule_length = 0.0
    seq_len = length(schema.sequence)
    # make sure the molecule is greater than 650
    while !(650 < molecule_length <= seq_len)
        molecule_length = rng.exponential(scale = schema.mol_length)
    end
    molecule_length = trunc(Int, molecule_length)
    adjusted_end = seq_len
    if !(schema.circular)
        # set the max start position to be (length - mol_length) to avoid overflow
        adjusted_end -= molecule_length
    end

    # if a singleton proportion is provided, conditionally drop the number of reads to 1
    if params.singletons > 0.0 && params.rng.uniform(0,1) <= params.singletons
        N = 1.0
    else
        if params.molecule_coverage < 1
            # set a minimum number of reads to 2 to avoid singletons
            _n = max(2.0, ((seq_len * params.molecule_coverage)) / (params.length_R1 + params.length_R2))
        else
            _n = params.molecule_coverage
        end
        # draw N from an exponential distribution with a minimum set to 2 reads to avoid singletons
        N = max(2.0, rng.exponential(_n))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, seq_len/(params.length_R1 + params.length_R2))
    end
    #TODO CREATE N inserts by sampling a normal dist or something like that
    
    #TODO STILL PYTHON CODE
    # this needs to iterate to a fixed number of times to try to create a molecule below a particular N
    # percentage else exit the program entirely
    for i in 1:attempts
        start_pos = params.rng.uniform(low = 0, high = adjusted_end)
        end_pos = start + molecule_length
        # 10 attempts to create N reads with <95% ambiguous bases
        for i in 1:10
       
        end
    end
end

function extract_sequence(sequence:::LongDNA{4}, breakpoints::Vector{UnitRange{Int}}, r1_len::Int, r2_len::Int)
    for breakpoint in breakpoints
        R1 = sequence[breakpoint.start:breakpoint.start+r1_len]
        R2 = sequence[breakpoint.stop:-1:breakpoint.stop-r2_len]
        if count(==(DNA_N), R1)/r1_len > 0.05 || count(==(DNA_N), R2)/r2_len > 0.05
            error("Too many Ns")
        end
    end

end

# Even more efficient version that maintains sorted intervals for faster overlap checking
function sample_nonoverlapping_ranges_fast(number_range, lengths)
    range_start = first(number_range)
    range_end = last(number_range)
    range_size = range_end - range_start + 1
    
    if sum(lengths) > range_size
        error("Cannot fit ranges of total length $(sum(lengths)) in range of size $range_size")
    end
    
    # Keep intervals sorted for efficient overlap checking
    occupied_intervals = Tuple{Int, Int}[]
    result_ranges = Vector{UnitRange{Int}}(undef, length(lengths))
    
    for (i, len) in enumerate(lengths)
        placed = false
        attempts = 0
        #TODO add attempts max
        max_attempts = min(1000, range_size)
        
        while !placed && attempts < max_attempts
            attempts += 1
            start_pos = rand(range_start:(range_end - len + 1))
            end_pos = start_pos + len - 1
            
            # Binary search for insertion point and overlap check
            insert_idx = searchsortedfirst(occupied_intervals, (start_pos, end_pos))
            
            # Check overlap with previous interval
            prev_overlaps = insert_idx > 1 && occupied_intervals[insert_idx-1][2] >= start_pos
            
            # Check overlap with next interval  
            next_overlaps = insert_idx <= length(occupied_intervals) && occupied_intervals[insert_idx][1] <= end_pos
            
            if !prev_overlaps && !next_overlaps
                # Insert the new interval in sorted order
                insert!(occupied_intervals, insert_idx, (start_pos, end_pos))
                result_ranges[i] = start_pos:end_pos
                placed = true
            end
        end
        
        if !placed
            error("Could not place range of length $len after $max_attempts attempts")
        end
    end
    
    return result_ranges
end

# Example usage comparing both methods:
number_range = 1:100
lengths = [10, 15, 8]



println("\n=== Fast Method (with sorted intervals) ===")
@time sampled_ranges2 = sample_nonoverlapping_ranges_fast(number_range, lengths)
println("Sampled ranges:")
for (i, r) in enumerate(sampled_ranges2)
    println("  Range $i: $r (length: $(length(r)))")
end

# Verification function
function check_overlaps(ranges)
    for i in 1:length(ranges)
        for j in (i+1):length(ranges)
            if !isempty(intersect(ranges[i], ranges[j]))
                return false
            end
        end
    end
    return true
end

println("\nVerification:")
println("Sequential method - No overlaps: $(check_overlaps(sampled_ranges1))")
println("Fast method - No overlaps: $(check_overlaps(sampled_ranges2))")

#=
def create_long_molecule(schema: Schema, rng, barcode: str, outputbarcode: str, wgsimparams, attempts: int) -> LongMoleculeRecipe|None:
    '''
    Randomly generates a long molecule and writes it to a FASTA file.
    Length of molecules is randomly distributed using an exponential distribution, with a minimum of 650bp.
    Returns a LongMoleculeRecipe that contains all the necessary information to simulate reads from that molecule.
    '''
    molnumber = getrandbits(32)
    len_interval = schema.end+1 - schema.start
    # make sure to cap the molecule length to the length of the interval/chromosome
    molecule_length = 0
    # make sure the molecule is greater than 650
    while molecule_length < 650 or molecule_length > len_interval:
        molecule_length = rng.exponential(scale = schema.mol_length)

    molecule_length = int(molecule_length)
    if schema.is_circular:
        # don't subtract the molecule length from the end, allow it to overflow since the sequence is repeated
        adjusted_end = len_interval
    else:
        # set the max start position to be (length - mol_length) to avoid overflow
        adjusted_end = len_interval - molecule_length
    # this needs to iterate to a fixed number of times to try to create a molecule below a particular N
    # percentage else exit the program entirely
    for i in range(attempts + 1):
        start = int(rng.uniform(low = 0, high = adjusted_end))
        end = start + molecule_length - 1
        fasta_seq = schema.sequence[start:end+1]
        N_count = fasta_seq.count('N')
        N_ratio = N_count/molecule_length
        if N_ratio < 0.7:
            normalized_length = molecule_length - N_count
            break
        elif i == attempts:
            return None

    if schema.is_circular:
        end %= len_interval
    # if a singleton proportion is provided, conditionally drop the number of reads to 1
    if schema.singletons > 0 and rng.uniform(0,1) <= schema.singletons:
        N = 1
    elif schema.mol_coverage < 1:
        # set a minimum number of reads to 2 to avoid singletons
        N = max(2, normalized_length * schema.mol_coverage)/(schema.read_length*2)
    else:
        # draw N from an exponential distribution with a minimum set to 2 reads to avoid singletons
        N = max(2, rng.exponential(schema.mol_coverage))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, normalized_length/(schema.read_length*2))
=#