"""
Takes a file with barcodes and validates them to be ATGC nucleotides and barcodes same length.
Returns nothing.
"""
function validate_barcodes(bc_list::Vector{LongSequence{DNAAlphabet{2}}})
    expected_length = length(first(bc_list))
    @inbounds for i in bc_list
        if length(i) != expected_length
            println("The barcodes provided must all be the same length " , i)
            break
        end
    end
    return
end


"""
Randomly generate `bc_count` barcodes, each with a length `bp` derived as
being the first number of bp to accommodate that many barcodes. Writes
to `filename`. Returns nothing.
"""
function generate_random_barcodes(bc_count::Int, segments::Int, filename::String)
    bp_max = [4^i for i in 4:30]
    n_bp = findfirst(x -> isless(bp_count, x), bp_max)
    _count = 0
    if segments == 1
        _gen = (string(shuffle(LongDNA{2}("ATCG"))) for i in 1:n_bp)
        generator = Iterators.Stateful(Iterators.product(_gen...))
    else
        #TODO is the # of barcodes really that important?
    end
        open(filename, "w") do file
        while _count < bc_count
            bc = reduce(*, iterate(generator)[1])
            println(file, bc)
            _count += 1
        end
    end
end

"""
Takes the barcode file and reads it line by line. Performs barcode validations and returns
a `BarcodeManifest`.
"""
function interpret_barcodes(infile::String, segments::Int, output_type::String)
    bc = Nothing
    try
        _file = safe_read(infile)
        try
            bc = unique(map(LongDNA{2}, readlines(_file)))
        catch
            println("The barcode file $infile contains characters that aren't valid nucleotides ATCG")
            return
        finally
            close(_file)
        end
    catch
        println("Cannot open $infile for reading because it\'s not recognized as either a plaintext or gzipped file.")
        return
    end
    validate_barcodes(bc)
    if segments == 1
        return BarcodeManifest(bc, output_type, length(bc))
    else
        return BarcodeManifest(Base.Iterators.product((bc for i in 1:segments)...), output_type, length(bc)^segments)
    end
end

function setup_barcode_output(output_type::String)
    if output_type == "haplotagging"
        #ints = [lpad(i,2,'0') for i in 1:96]
        return Base.Iterators.Stateful(
            Base.Iterators.product('A', 1:96, 'C', 1:96, 'B', 1:96, 'D', 1:96)
        )
    elseif output_type == "stlfr"
        return Base.Iterators.Stateful(
            Base.Iterators.product(1:2000, 1:2000, 1:2000)
        )
    else
        return Nothing
    end
end

"""
Format a haplotagging style barcode for output
"""
function format_output_barcode(bc::Tuple{Tuple{Char, Int64, Char, Int64, Char, Int64, Char, Int64}, Nothing}, standard::Bool = true)::String
    formatted = "BX:Z:" * mapreduce(*, bc[1]) do x
        if x isa Int
            lpad(x, 2, "0")
        else
            x
        end
    end
    standard ? "$formatted\tVX:i:1" : formatted
end

"""
Format a stLFR style barcode for output
"""
function format_output_barcode(bc::Tuple{Tuple{Int64, Int64, Int64}, Nothing}, standard::Bool = true)::String
    formatted = join(bc[1],"_")
    standard ? "BX:Z:$formatted\tVX:i:1" : "#$formatted" 
end

"""
Format a TELLseq/10X style barcode for output (i.e. just return it)
"""
function format_output_barcode(bc::String, standard::Bool = true)::String
    standard ? "BX:Z:$bc\tVX:i:1" :  bc
end

"""
Format a nucleotide barcode created by iterating the generator
"""
function format_nuc_barcode(bc)::String
    join(bc[1],"")
end