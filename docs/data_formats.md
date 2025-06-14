# Linked-read Data Formats

Linked-read chemistries continue to evolve and it seems like every new method wants to use their own bespoke
convention for putting barcodes in FASTQ files (because that's exactly what's happening, unfortunately). Until
such a day comes when linked-read FASTQ data formats will be properly standardized, we are sort of forced to just let
tyrants roam about with their hubris. It is for this very frustrating reason that Mimick supports different
input and output linked-read types.

You can specify the linked-read barcode chemistry to simulate via `--lr-type` as well as
the output format of FASTQ files (default is the same as barcode type). For example, you
can generate 96 barcodes (common haplotagging style), select `--lr-type stlfr`
(combinatorial 3-barcode on R2 read), and have `--output-type tellseq` (`@seqid:barcode` header format).

## Linked-read simulation types
These are the options available to `--lr-type`. Each has a different approach to using the
supplied/generated barcode list and the length of the resulting reads to compensate for where
the barcode(s) would have been attached to.

| --lr-type       | Format                                                 |
|:----------------|:-------------------------------------------------------|
| `10x`/`tellseq` | single barcode on R1                                   |
| `haplotagging`  | R1 and R2 each have different combinatorial 2-barcodes |
| `stlfr`         | combinatorial 3-barcode on R2                          |

In practice, this means that if you input 400,000  16bp barcodes for `10x` or `tellseq` simulation, the
resulting R1 reads will be 16bp shorter, because that's where that barcodes would have been. 
Unlike `10x` and `tellseq`, which use barcodes directly, you need far fewer barcodes as input for
`haplotagging` and `stlfr`. For example, standard haplotagging uses 96 barcodes per segment and standard
stlfr uses 1537 barcodes per segment. Haplotagging will make $N^4$ barcode combinations, whereas stLFR
will make $N^3$ combinations.

## Linked-read output types
Like discussed above, there are _options_ for how the resulting linked-read data can look. Why would you want one
format over another? Well, it could be personal preference or the software you want to use is configured for a very
specific format (which is a **problem** for the linked-read ecosystem). Regardless of the _kind_ of linked-read
experiment you chose with `--lr-type`, you can specify any of the linked-read types as the output format with `--output-type`.

| --output-type  | Barcode Location                                      | Example                    |
|:---------------|:------------------------------------------------------|:---------------------------|
| `10x`          | start of R1 sequence                                  | `ATAGACCATAGA`GGACA...     |
| `haplotagging` | sequence header as `BX:Z:ACBD`                        | `@SEQID BX:Z:A0C331B34D87` |
| `standard`     | sequence header as `BX:Z:BARCODE`, no specific format | `@SEQID BX:Z:ATACGAGACA`   |
| `stlfr`        | appended to sequence ID via `#1_2_3`                  | `@SEQID#1_354_39`          |
| `tellseq`      | appended to sequence ID via `:ATCG`                   | `@SEQID:TATTAGCAC`         |

## Proper read pairing
A consideration for downstream application using the simulated linked reads is that you might need to
make sure the paired-end reads are properly paired, since Mimick does not enforce validating proper
read pairs from `wgsim`. Sometimes software gets fussy when it expects proper (synchronized) reads pairs
and doesn't get them, like `bwa`, for example. A quick way to accomplish this is using `seqkit` (or similar):

```bash
# if the forward/reverse specification is the /1 /2 format
seqkit pair --id-regexp '^(\S+)\/[12]' -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz

# if the forward/reverse specification is the modern 1:N:0:ATTACA format
seqkit pair -1 sample_0${i}.R1.fq.gz -2 sample_0${i}.R2.fq.gz
```
You can optionally use the `-u` flag to ask `seqkit` to also save the unpaired reads.