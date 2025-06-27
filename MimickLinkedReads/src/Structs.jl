
struct SimParams
    error::Float32
    #mutation::Float32
    #indels::Float32
    #extindels::Float32
    read_distance::Int32
    distance_stdev::Float32
    length_R1::Int16
    length_R2::Int16
    randomseed::Int64
    output_dir::String
    prefix::String
    function SimParams(error, read_distance, distance_stdev, length_R1, length_R2, randomseed, prefix)
        _prefix = Base.Filesystem.basename(prefix)
        _outdir = Base.Filesystem.dirname(prefix)
        return new(error, read_distance, distance_stdev, length_R1, length_R2, randomseed, _outdir, _prefix)
    end
end

function Base.show(io::IO, data::SimParams)
    println("Object of type SimParams")
    for i in fieldnames(SimParams)
        val = getfield(data, i)
        println(io, " $i::", typeof(val), " ", val)
    end
end




