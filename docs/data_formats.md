# Linked-read Data Formats

Linked-read chemistries continue to evolve and it seems like every new method wants to use their own bespoke
convention for putting barcodes in FASTQ files (because that's exactly what's happening, unfortunately). Until
such a day comes when linked-read FASTQ data formats will be properly standardized, we are sort of forced to just let
tyrants roam about with their hubris. It is for this very frustrating reason that Mimick supports different
input and output linked-read types.

You can specify the linked-read barcode chemistry to simulate via `--lr-type` as well as
the output format of FASTQ files (default is the same as barcode type). For example, you
can generate 96 barcodes (common haplotagging style), select `--barcode-type stlfr`
(combinatorial 3-barcode on R2 read), and have `--output-format tellseq` (`@seqid:barcode` header format).

## Linked-read simulation types
These are the options available to `--lr-type`. Each has a different approach to using the
supplied/generated barcode list and the length of the resulting reads to compensate for where
the barcode(s) would have been attached to.

| --lr-type | Format |
|:------------------|:-------|
|`10x`/`tellseq`   | single barcode on R1 |
|`haplotagging`  | R1 and R2 each have different combinatorial 2-barcodes |
|`stlfr`         | combinatorial 3-barcode on R2 |

In practice, this means that if you input 400,000  16bp barcodes for `10x` or `tellseq` simulation, the
resulting R1 reads will be 16bp shorter, because that's where that barcodes would have been. 
Unlike `10x` and `tellseq`, which use barcodes directly, you need far fewer barcodes as input for
`haplotagging` and `stlfr`. For example, standard haplotagging uses 96 barcodes per segment and standard
stlfr uses 1537 barcodes per segment. Haplotagging will make $N^4$ barcode combinations, whereas stLFR
will make $N^3$ combinations. As a result, you will need far fewer barcodes as input for these two methods.

!> Haplotagging typically uses 96 barcodes. stLFR typically uses 1537 barcodes. Using too many barcodes for the combinatorial methods may result in an error

## Linked-read output types
Like discussed above, there are _options_ for how the resulting linked-read data can look. Why would you want one
format over another? Well, it could be personal preference or the software you want to use is configured for a very
specific format (which is a **problem** for the linked-read ecosystem). Regardless of the _kind_ of linked-read
experiment you chose with `--lr-type`, you can specify any of the linked-read types as the output format with `--output-format`.

| --output-type | Barcode Location | Example |
|:-----------------|:-------|:---------------------|
|`10x`           | start of R1 sequence | `ATAGACCATAGA`GGACA... |
|`haplotagging`  | sequence header as `BX:Z:ACBD` |  `@SEQID BX:Z:A0C331B34D87` |
|`standard`      | sequence header as `BX:Z:BARCODE`, no specific format | `@SEQID BX:Z:ATACGAGACA` |
|`stlfr`         | appended to sequence ID via `#1_2_3` | `@SEQID#1_354_39` |
|`tellseq`       | appended to sequence ID via `:ATCG` | `@SEQID:TATTAGCAC` |