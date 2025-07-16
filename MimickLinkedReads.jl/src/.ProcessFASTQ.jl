
function start_writer_thread(writer::FastqWriter, format::Symbol)
    return @spawn begin
        FASTQ.Writer(GzipCompressorStream(open(writer.prefix * ".R1.fq.gz", "w"))) do R1; FASTQ.Writer(GzipCompressorStream(open(writer.prefix * ".R2.fq.gz", "w"))) do R2 
            while writer.running || isready(writer.queue)
                try
                    # Wait for items with timeout to allow checking running status
                    if isready(writer.queue)
                        molecule = take!(writer.queue)
                        if molecule === :shutdown
                            break
                        end
                        # convert to fastq records and write R1
                        for record in format_R1(format, molecule)
                            write(R1, record)
                        end
                        flush(R1)

                        # convert to fastq records and write R2
                        for record in format_R2(format, molecule)
                            write(R2, record)
                        end
                        flush(R2)
                    else
                        # Small sleep to prevent busy waiting
                        sleep(0.00001)
                    end
                catch e
                    if isa(e, InvalidStateException) && !writer.running
                        # Queue was closed, exit gracefully
                        break
                    else
                        println(stderr, "Error in writer thread: ", e)
                    end
                end
            end
        end; end
        println("Writer thread stopped")
    end
end


"""
Submit a ProcessedMolecule to the writing thread to be converted to FASTQ records and written
"""
function submit!(writer::FastqWriter, molecule::ProcessedMolecule)
    if writer.running
        put!(writer.queue, molecule)
    else
        throw(ArgumentError("Writer is not running"))
    end
end

"""
Gracefully stop the writer thread
"""
function stop!(writer::FastqWriter)
    if writer.running
        writer.running = false
        # Signal to stop
        put!(writer.queue, :shutdown)
        # Wait for thread to finish
        wait(writer.writer_task)
        close(writer.queue)
    end
end

"""
Check if writer is still running
"""
function is_running(writer::FastqWriter)
    return writer.running && !istaskdone(writer.writer_task)
end
