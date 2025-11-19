---
label: Linked-Read Data Formats
icon: file-binary
---

Linked-read chemistries continue to evolve and it seems like every new method wants to use their own bespoke
convention for putting barcodes in FASTQ files (because that's exactly what's happening, unfortunately). Until
such a day comes when linked-read FASTQ data formats will be properly standardized, we are sort of forced to just let
tyrants roam about with their hubris. It is for this very frustrating reason that Mimick supports different
input and output linked-read types.

## Linked-read simulation types
You can specify the linked-read barcode chemistry to simulate using the combination of `--read-lengths`, and
`--format`. For example, you can simulate the common stLFR style (combinatorial 3-barcode on R2)
with `--format stlfr` and `--read-lengths 150,108`. You can also mix-match these options, such as
`--format haplotagging --read-lengths 134,150` (`@seqid:barcode` header format).

The table below serves as a guide for the configurations for the common linked-read varieties: 

{.compact}
| Chemistry    | `--read-lengths` | Format                                  | barcode `--format`      |
|:-------------|:-----------:|:---------------------------------------------|:------------------------|
| 10x          |  `134,150`  | single barcode on R1                         | `tellseq`               |
| tellseq      |  `132,150`  | single barcode on R1                         | `tellseq`               |
| haplotagging |  `150,150`  | I1 and I2 each with combinatorial 2-barcodes | `standard:haplotagging` |
| stlfr        |  `150,108`  | combinatorial 3-barcode on R2                | `stlfr`                 |

The simulation process never actually includes the barcodes in the reads, so the read lengths you
specify with `--lengths` will be the **final demultiplexed read lengths**. Unlike 10x and tellseq, which use barcodes directly,
far fewer unique barcodes are needed for the combinatorial chemistries (haplotagging and stlfr). For example, standard haplotagging uses 96
barcodes per segment and standard stlfr uses 1537 barcodes per segment. Haplotagging will make $N^4$ barcode combinations,
whereas stLFR will make $N^3$ combinations.

## Linked-read output types
Like discussed above, there are _options_ for how the resulting linked-read data can look. Why would you want one
format over another? Well, it could be personal preference or the software you want to use is configured for a very
specific format (which is a **problem** for the linked-read ecosystem). Regardless of the _kind_ of linked-read
experiment you are trying to do, you can specify any of the linked-read types as the output format with `--format`.
You can suffix `standard` with `:haplotagging` or `:stlfr` (e.g. `standard:stlfr`) to output the standard format
with that kind of barcode encoding style, otherwise `standard` (no suffix) will use the nucleotide barcode.

{.compact}
| --format    | Barcode Location                         | Example                    |
|:-----------------|:-----------------------------------------|:---------------------------|
| `10x`            | start of R1 sequence                     | `ATAGACCATAGA`GGACA...     |
| `haplotagging`   | sequence header as `BX:Z:ACBD`           | `@SEQID BX:Z:A0C331B34D87` |
| `standard[:...]` | sequence header as `BX:Z:BARCODE VX:i:N` | `@SEQID BX:Z:ATACGAGACA`   |
| `stlfr`          | appended to sequence ID via `#1_2_3`     | `@SEQID#1_354_39`          |
| `tellseq`        | appended to sequence ID via `:ATCG`      | `@SEQID:TATTAGCAC`         |

## Proper read pairing
Mimick should properly pair reads, but in the event you need to do that manually, you can use `seqkit`:
```bash
# if the forward/reverse specification is the /1 /2 format
seqkit pair --id-regexp '^(\S+)\/[12]' -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz

# if the forward/reverse specification is the modern 1:N:0:ATTACA format
seqkit pair -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz
```
You can optionally use the `-u` flag to ask `seqkit` to also save the unpaired reads.