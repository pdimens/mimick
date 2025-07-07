
"""
Given a `ProcessedMolecule`, will convert the sequences therein into formatted FASTQ records and write them to
appropriate R1 and R2 files. 
"""
function convert2fastq(molecule::ProcessedMolecule, qual::String)
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    # R1 reads
    #TODO LOGIC FOR FORWARD/REVERSE NOTATION
    #TODO LOGIC FOR BARCODE FORMAT
    FASTQ.Writer(open("some_file.R1.fq", "w")) do R1
        @inbounds _qual = qual^length(molecule.read_sequences.first[begin])
        @inbounds for (breakpoint,seq) in zip(molecule.read_breakpoints, molecule.read_sequences.first)
            write(R1, FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/1", seq, _qual))
        end
    end
    FASTQ.Writer(open("some_file.R2.fq", "w")) do R2
        @inbounds _qual = qual^length(molecule.read_sequences.second[begin])
        @inbounds for (breakpoint,seq) in zip(molecule.read_breakpoints, molecule.read_sequences.second)
            write(R2, FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/2", seq, _qual))
        end
    end
end