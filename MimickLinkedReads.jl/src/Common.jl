


function safe_read(filename::String)
    try
        readline(GzipDecompressorStream(open(filename, "r")))
        return GzipDecompressorStream(open(filename, "r"))
    catch e
        if isa(e, CodecZlib.ZlibError) || isa(e, EOFError)
            return open(filename, "r")
        else
            rethrow(e)
        end
    end
end