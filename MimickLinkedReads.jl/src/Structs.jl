
struct SimParams
    output_dir::String
    prefix::String
    error::Float32
    #mutation::Float32
    #indels::Float32
    #extindels::Float32
    read_distance::Int32
    distance_stdev::Float32
    length_R1::Int16
    length_R2::Int16
    molecule_length::Int
    molecule_coverage::Float64
    singletons::Float64
    circular::Bool
    #TODO rng
    randomseed::Int64
    function SimParams(prefix, error, read_distance, distance_stdev, length_R1, length_R2, molecule_length, molecule_coverage, singletons, circular; randomseed::Int)
        _prefix = Base.Filesystem.basename(prefix)
        _outdir = Base.Filesystem.dirname(prefix)
        return new(_outdir, _prefix, error, read_distance, distance_stdev, length_R1, length_R2, molecule_length, molecule_coverage, singletons, circular, randomseed)
    end
end

function Base.show(io::IO, data::SimParams)
    println("Object of type SimParams")
    for i in fieldnames(SimParams)
        val = getfield(data, i)
        println(io, " $i::", typeof(val), " ", val)
    end
end




