# Linked-read Data Formats

Linked-read chemistries continue to evolve and it seems like every new method wants to use their own bespoke
convention for putting barcodes in FASTQ files (because that's exactly what's happening, unfortunately). Until
such a day comes when linked-read FASTQ data formats will be properly standardized, we are sort of forced to just let
tyrants roam about with their hubris. It is for this very frustrating reason that Mimick supports different
input and output linked-read types.

## Linked-read simulation types
You can specify the linked-read barcode chemistry to simulate using the combination of `--lengths`, `--segments` and
`--output-type`. For example, you can generate 96 6bp barcodes (common haplotagging style), select `--segments 3`
(combinatorial 3-barcode, the stLFR style), and have `--output-type tellseq` (`@seqid:barcode` header format).

The table below serves as a guide for the configurations for the common linked-read varieties: 

| Chemistry    | `--segments` | `--lengths` | Format                                       | `--output-type` default |
|:-------------|:------------:|:-----------:|:---------------------------------------------|:------------------------|
| 10x/tellseq  |     `1`      |  `134,150`  | single barcode on R1                         | `tellseq`               |
| haplotagging |     `4`      |  `150,150`  | I1 and I2 each with combinatorial 2-barcodes | `standard:haplotagging` |
| stlfr        |     `3`      |  `150,108`  | combinatorial 3-barcode on R2                | `stlfr`                 |

The simulation process never actually includes the barcodes in the reads, so the read lengths you
specify with `--lengths` will be the **final demultiplexed read lengths** (this means `--output-type 10x` will
have a longer R1 to add the barcode inline). Unlike 10x and tellseq, which use barcodes directly,
you need far fewer barcodes as input for haplotagging and stlfr. For example, standard haplotagging uses 96
barcodes per segment and standard stlfr uses 1537 barcodes per segment. Haplotagging will make $N^4$ barcode combinations,
whereas stLFR will make $N^3$ combinations.

## Linked-read output types
Like discussed above, there are _options_ for how the resulting linked-read data can look. Why would you want one
format over another? Well, it could be personal preference or the software you want to use is configured for a very
specific format (which is a **problem** for the linked-read ecosystem). Regardless of the _kind_ of linked-read
experiment you are trying to do, you can specify any of the linked-read types as the output format with `--output-type`.
You can suffix `standard` with `:haplotagging` or `:stlfr` (e.g. `standard:stlfr`) to output the standard format
with that kind of barcode encoding style, otherwise `standard` (no suffix) will use the nucleotide barcode.

| --output-type    | Barcode Location                         | default for | Example                    |
|:-----------------|:-----------------------------------------|:-----------:|:---------------------------|
| `10x`            | start of R1 sequence                     |             | `ATAGACCATAGA`GGACA...     |
| `haplotagging`   | sequence header as `BX:Z:ACBD`           |             | `@SEQID BX:Z:A0C331B34D87` |
| `standard[:...]` | sequence header as `BX:Z:BARCODE VX:i:N` | all others  | `@SEQID BX:Z:ATACGAGACA`   |
| `stlfr`          | appended to sequence ID via `#1_2_3`     | `-x 3`      | `@SEQID#1_354_39`          |
| `tellseq`        | appended to sequence ID via `:ATCG`      | `-x 1`      | `@SEQID:TATTAGCAC`         |

## Proper read pairing
A consideration for downstream application using the simulated linked reads is that you might need to
make sure the paired-end reads are properly paired, since Mimick does not enforce validating proper
read pairs from `wgsim`. Sometimes software gets fussy when it expects proper (synchronized) reads pairs
and doesn't get them, like `bwa`, for example. Often times the reads are already properly paired and this isn't
a concern. If you need to ensure proper pairing, a quick way to do it is using `seqkit` (or similar):

```bash
# if the forward/reverse specification is the /1 /2 format
seqkit pair --id-regexp '^(\S+)\/[12]' -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz

# if the forward/reverse specification is the modern 1:N:0:ATTACA format
seqkit pair -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz
```
You can optionally use the `-u` flag to ask `seqkit` to also save the unpaired reads.