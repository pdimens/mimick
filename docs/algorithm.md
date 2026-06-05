---
label: Algorithm Overview
icon: git-compare
---

# The simulation process
Or _algorithm_, if you want to be specific. Note that mentions of "contig" below are synonymous with
"chromosome" if your input FASTA files are chromosome-scale.

Mimick 2.0 had a complete overhaul that rewrote the entire simulation process. Then the Mimick 3.0 Julia
port rewrote the algorithm again, making it *much* more streamlined. Still missing is sequence errors modeled after real Illumina data. That's on the wish list.

## Single-Sample

### 0: Parse the FASTA file(s)
The FASTA files are processed into schema that contain characteristics for each contig, along with the
sequence in it. These schema are part of the recipes that get processed during simulation. If `--mask-n` was
used, the `sequence` will have instances of `N` replaced with a random `ATCG` base. Schemas look like this:

```julia
struct Schema
    haplotype::Int
    chrom::String
    tracker::SchemaTracker
    sequence::LongSequence{DNAAlphabet{4}}
```
If you aren't familiar with Julia code, the syntax is `name::Type`. The `SchemaTracker` is just two integers tracking the goal and the current total:
```julia
mutable struct SchemaTracker
    reads_current::Int
    reads_required::Int
```

What follows is all done in a loop
### 1: Pick a barcode
Barcodes get created using a generator depending on the `--format`.
You can think of the barcode as the oligo-covered bead in TELLseq/stLFR/haplotagging
chemistries.

```python
# examples
haplotagging = "A43C21B01D93"
stlfr = "1_32_441"
tellseq = "ATAGGGACACAGATGACC"
10x = "TAGCGAACAGATGACC"
```

### 2: Making the molecule recipes
Based on the `--molecules-per` parameter, Mimick will decide how many unrelated molecules will be associated with the barcode.
Then, it's a matter of randomly choosing (with replacement) which schema to make molecules from. Remember, each contig for each
haplotype has its own schema, so unrelated molecules can come from anywhere, which is a very important feature of these simulations.
The "molecules" are just recipes describing which contig/haplotype, which section of it, how many paired
reads need to be generated from it to meet `--molecule-coverage` targets, etc. There are a few criteria that need to be met for a molecule:
- must be at least 650bp
- if larger than the source contig, caps to its maximum possible length
    - e.g. a simulated molecule length was 2000, but the contig is only 1000bp, molecule becomes 1000bp 
- molecule end position cannot exceed that of the source contig

The number of reads are calculated based on `--molecule-coverage` and the molecule size, with a minimum of 2 to prevent singletons. If
a value >0 is set for `--singletons`, then there is that percent chance the number of paired reads are overwritten to be 1
(e.g. `--singletons 0.5` = 50% chance of singletons for each molecule). Molecules look like this:
```julia
mutable struct ProcessedMolecule
    haplotype::Int     # 1
    barcode::String    # "4_5_11"
    chrom::String      # "contig1"
    position::UnitRange{Int}   # 30000:45000
    read_breakpoints::Vector{UnitRange{Int}}   # [30090:30490, 43700:44000]
    read_sequences::Pair{Vector{String}, Vector{String}} # [[R1s] => [R2s]]
```

### 3: Creating sequences from a molecule
Once we have a "molecule" made and the number of reads to generate from it, an algorithm creates `n` random non-overlapping
inserts (start + end positions) along the molecule. Once the insert breakpoints are found, their sequences are finally extracted
from the haplotype at those breakpoints, and if they contain >5% `N` sequences, new insert positions are attempted. If `N` bases were masked with `--mask-n`, then there is no concern of too many `N`s.

### 4: Write FASTQ records
After an adequate set of sequences is determined from the molecule, the sequences get converted into proper FASTQ records in the
`--format` of your choosing. This is done using a BGZip stream for each R1/R2 file. When supplying `--threads` to single-sample
simulation, the threads are applied to the BGZip compressors because compression + writing to disk is the performance bottleneck.

!!-question Careful parallelization
It might seem like parallelizing read simulation would make more sense than file IO. I thought so too and that was the first implementation in the Julia code. It turned out that the code was so efficient that in many small-to-medium cases, the simulations were 100% finished before the first FASTQ record was written to disk and most of the run time was just waiting
for the compressed FASTQ records to be written to disk. Providing additional threads to the compression and write turned out
to actually be the faster solution.
!!!


### 5: Monitor schema targets
The number of reads that were generated for every molecule are tracked in the schema for each contig-haplotype. Once a schema reaches its reads target, that schema will be removed from the list of schema that are randomly sampled
to determine molecules in [Step 2](#2-making-the-molecule-recipes). This ensures that genomic `--coverage` is honored.

### 5: Repeat until all schema read targets are met
Once there are no more schema left to sample, simulation is done!

## Multi-Sample
This alternative method of using Mimick uses the same simulation algorithm as the single-sample approach, but is intended for multi-sample simulation, saving a lot of effort. The inputs are a single 
FASTA reference and a VCF file (can be vcf.gz) that contains _samples_ and SNPs/indels. The main difference is that this method
parallelizes over multiple samples. It does so by using the variant information in the VCF file to "paint" that sample's variants onto the FASTA reference and then proceeding as normal. In psuedo-code, that looks like:

```python
# pseudo-code
samples = get_samples(vcf)
for sample in samples:     <- parallelized
    schemas = []
    variants = get_variants(sample, vcf)
    for hap in 1:ploidy(variants):
        haplotype = paint_haplotype(variants, hap, fasta)
        append(schemas, haplotype)
    simulate_reads(schemas)
```