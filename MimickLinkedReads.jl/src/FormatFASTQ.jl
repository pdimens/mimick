## R1 ##
"""
`format_R1(format::Symbol, molecule::ProcessedMolecule)`

Converts a ProcessedMolecule into a Vector of FASTQRecords, formatted in the style of `format`
(i.e. `haplotagging`, `stlfr`, `tellseq`, `tenx`, `standard`).
"""
format_R1(format::Symbol, molecule::ProcessedMolecule) = format_R1(Val(format), molecule)

function format_R1(::Val{:stlfr}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.first[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.first)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)#$(molecule.barcode) 1:N:0:ATGACA", seq, _qual)
    end
end

function format_R1(::Val{:haplotagging}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.first[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.first)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/1\tBX:Z:$(molecule.barcode)", seq, _qual)
    end
end

function format_R1(::Val{:tellseq}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.first[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.first)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop):$(molecule.barcode)  1:N:0:ATGACA", seq, _qual)
    end
end

function format_R1(::Val{:tenx}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^(length(molecule.read_sequences.first[1]) + length(molecule.barcode))
    map(zip(molecule.read_breakpoints, molecule.read_sequences.first)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop) 1:N:0:ATGACA", molecule.barcode * seq, _qual)
    end
end

function format_R1(::Val{:standard}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.first[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.first)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/1\tVX:i:1\tBX:Z:$(molecule.barcode)", seq, _qual)
    end
end

## R2 ##
"""
`format_R2(format::Symbol, molecule::ProcessedMolecule)`

Converts a ProcessedMolecule into a Vector of FASTQRecords, formatted in the style of `format`
(i.e. `haplotagging`, `stlfr`, `tellseq`, `tenx`, `standard`).
"""
format_R2(format::Symbol, molecule::ProcessedMolecule) = format_R2(Val(format), molecule)

function format_R2(::Val{:stlfr}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.second[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.second)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)#$(molecule.barcode) 2:N:0:ATGACA", seq, _qual)
    end
end

function format_R2(::Val{:haplotagging}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.second[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.second)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/2\tBX:Z:$(molecule.barcode)", seq, _qual)
    end
end

function format_R2(::Val{:tellseq}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.second[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.second)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop):$(molecule.barcode)  2:N:0:ATGACA", seq, _qual)
    end
end

function format_R2(::Val{:tenx}, molecule::ProcessedMolecule)::Vector{FASTQRecord}
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.second[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.second)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop) 2:N:0:ATGACA", seq, _qual)
    end
end

function format_R2(::Val{:standard}, molecule::ProcessedMolecule)
    header = "$(molecule.chrom):$(molecule.haplotype)|$(molecule.position.start):$(molecule.position.stop)|"
    @inbounds _qual = "I"^length(molecule.read_sequences.second[1])
    map(zip(molecule.read_breakpoints, molecule.read_sequences.second)) do (breakpoint,seq)
        FASTQRecord(header * "$(breakpoint.start):$(breakpoint.stop)/2\tVX:i:1\tBX:Z:$(molecule.barcode)", seq, _qual)
    end
end