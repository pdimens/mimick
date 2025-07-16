## R1 ##
"""
`format_R1(format::Symbol, molecule::ProcessedMolecule)`

Converts a ProcessedMolecule into a String of FASTQ records, formatted in the style of `format`
(i.e. `haplotagging`, `stlfr`, `tellseq`, `tenx`, `standard`).
"""
format_R1(format::Symbol, molecule::ProcessedMolecule)::String = format_R1(Val(format), molecule)

function format_R1(::Val{:stlfr}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.first[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])#" * molecule.barcode * " 1:N:0:ATGACA\n" * molecule.read_sequences.first[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R1(::Val{:haplotagging}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.first[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])/1\tBX:Z:" * molecule.barcode * "\n" * molecule.read_sequences.first[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R1(::Val{:tellseq}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.first[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i]):" * molecule.barcode * " 1:N:0:ATGACA\n" * molecule.read_sequences.first[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R1(::Val{:tenx}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^(length(molecule.read_sequences.first[1]) + length(molecule.barcode))
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i]) 1:N:0:ATGACA\n" * molecule.barcode * molecule.read_sequences.first[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R1(::Val{:standard}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.first[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])/1\tVX:i:1\tBX:Z:" * molecule.barcode * "\n" * molecule.read_sequences.first[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

## R2 ##
"""
`format_R2(format::Symbol, molecule::ProcessedMolecule)`

Converts a ProcessedMolecule into a String of FASTQ records, formatted in the style of `format`
(i.e. `haplotagging`, `stlfr`, `tellseq`, `tenx`, `standard`).
"""
format_R2(format::Symbol, molecule::ProcessedMolecule)::String = format_R2(Val(format), molecule)

function format_R2(::Val{:stlfr}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.second[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])#" * molecule.barcode * " 2:N:0:ATGACA\n" * molecule.read_sequences.second[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R2(::Val{:haplotagging}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.second[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])/2\tBX:Z:" * molecule.barcode * "\n" * molecule.read_sequences.second[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R2(::Val{:tellseq}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.second[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i]):" * molecule.barcode * " 2:N:0:ATGACA\n" * molecule.read_sequences.second[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R2(::Val{:tenx}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.second[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i]) 2:N:0:ATGACA\n" * molecule.read_sequences.second[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end

function format_R2(::Val{:standard}, molecule::ProcessedMolecule)::String
    @inbounds begin
        header = "@" * molecule.chrom * ":$(molecule.haplotype)|$(molecule.position)|"
        fq = ""
        _qual = "I"^length(molecule.read_sequences.second[1])
        for i in eachindex(molecule.read_breakpoints)
            fq *= header * "$(molecule.read_breakpoints[i])/2\tVX:i:1\tBX:Z:" * molecule.barcode * "\n" * molecule.read_sequences.second[i] * "\n+\n" * _qual * "\n"
        end
    end
    return fq
end