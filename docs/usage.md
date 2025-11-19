---
label: Using Mimick
icon: terminal
order: 98
---

```bash usage
mimick options... BARCODES FASTA1 FAST2...
```
Use `--help` or `-h` or call `mimick` without arguments to call up the docstring.

The minimum input files required by Mimick is a single FASTA file, uncompressed or gzip compressed.
The result will be a single set of paired-end reads, a GFF file of mutations, and a manifest of all the molecules created.
Mimick scales with the number of threads provided, although not linearly (the slowest part is writing compressed files to disk).

If you want to separate out the haplotypes from the final output, you can
leverage the fact that all the simulated reads start with `@HAP:X_` where `X` is the haplotype number, starting with `1`,
corresponding to the order in which FASTA files were provided. The read names also
include source contig names and other identifying features that can be extracted similarly.

```bash
for i in {1..2}; do
    zgrep -A3 \"^@HAP:X_\" output_prefix.R$i.fq.gz | gzip > out.R$i.fq.gz
done
```


==- All Options
For completeness, the table below is all the command line arguments and options

{.compact}
| short name | long name             |                         default                          | description                                                                                                                                    | notes                                                                               |
|:----------:|:----------------------|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
|            | `FASTA`               |                                                          | input fasta file(s)                                                                                                                            | REQUIRED                                                                            |
|    `-c`    | `--circular`          |                                                          | toggle to let Mimick know the input contigs are circular/prokaryotic                                                                           |                                                  |
|    `-o`    | `--output-prefix`     |                     `simulated/`                         | output file prefix                                                                                                                             |                                                                                     |
|    `-f`    | `--format`            | [varies](./data_formats.md#linked-read-simulation-types) | output format of FASTQ files                                                                                                                   |                                                                                     |
|    `-q`    | `--quiet`             |                                                          | toggle to hide the progress bar                                                                                             |                                                                                     |
|    `-s`    | `--seed`              |                                                          | random seed for simulation                                                                                                                     | optional and useful for reproducibility                                             |
|    `-t`    | `--threads`           |                           `2`                            | number of threads to use for simulation                                                                                                        |                                                                                     |
|    `-g`    | `--genomic-coverage`  |                           `30`                           | mean coverage target for simulated data                                                                                                        |                                                                                     |
|    `-i`    | `--insert-size`       |                          `500`                           | outer distance between the two ends in bp                                                                                                      | must be >`--read-lengths`                                                                 |
|    `-d`    | `--insert-stdev`      |                           `50`                           | standard deviation of --distance                                                                                                               |                                                                                     |
|    `-A`    | `--molecule-attempts` |                          `300`                           | number of attempts to create a molecule with <70% ambiguous bases before exiting with an error                                                 |                                                                                     |
|    `-C`    | `--molecule-coverage` |                          `0.2`                           | mean percent coverage per molecule if <1, else mean number of reads per molecule                                                               |                                                                                     |
|    `-L`    | `--molecule-length`   |                         `80000`                          | mean length of molecules in bp, drawn from exponential distribution                                                                            |                                                                                     |
|    `-N`    | `--molecules-per`     |                           `3`                            | mean number of unrelated molecules per barcode, drawn from an exponential distribution. If negative, (e.g. `-2`) will be fixed for that number |                                                                                     |
|    `-l`    | `--read-lengths`      |                        `150,150`                         | length of R1,R2 reads in bp                                                                                                                    | each must be >10 and separated by a comma, no spaces                                |
|    `-S`    | `--singletons`        |                           `0`                            | proportion of barcodes will only have a single read pair                                                                                       |
|    `-v`    | `--vcf`               |                                                          | VCF-formatted file containing genotypes from which to create per-sample haplotypes                                                             |                                              |

===

## Required Arguments
Running Mimick will always require a barcode specification and at least one input FASTA file.

+++ FASTA file(s)
You will need at least 1 fasta file as input, which goes at the very end of the command. If providing more than 1
FASTA for a non-haploid species, it's assumed that each fasta file is a different haplotype **of the same genome**.
There is no strict enforcement of the fasta files being from the same genome, it's just how you'll probably want to use
the simulator (but I'm open to being surprised). 

```bash fasta inputs
mimick hap1.fasta hap2.fasta
#      ^fasta 1   ^fasta 2
```
### Circular DNA
The `--circular` toggle/flag tells Mimick to treat each contig
within each FASTA as circular when creating molecules-- this is probably how you'll want to simulate prokaryotic/microbial
genomes. When using `--circular`, you'll see in the `.molecules` file that some molecule end positions may occur before the start
positions, which is an artifact of the circularization. Those positions reflect the start/end positions on the linear sequence in
the FASTA file and should be interpreted as "started at `start` and reached the end of the contig, then wrapped around to the
beginning of the contig and kept going until `end`." The corresponding length of the molecule will be accurate and the math
should make sense as well: $end = start + molecule\_size - contig\_size$

+++

## Options
+++ General options
These options control inputs/outputs and resources

{.compact}
| short name | long name         | default                                                  | description                                                          |
|:----------:|:------------------|:--------------------------------------------------------:|:---------------------------------------------------------------------|
|    `-c`    | `--circular`      | `False`                                                  | toggle to let Mimick know the input contigs are circular/prokaryotic |
|    `-o`    | `--output-prefix` | `simulated/`                                          | output file prefix                                                   |
|    `-f`    | `--format`        | [varies](./data_formats.md#linked-read-simulation-types) | output format of FASTQ files                                         |
|    `-q`    | `--quiet`         | `False`                                                  | toggle to hide progress bar                                          |
|    `-S`    | `--seed`          |                                                          | random seed for simulation                                           |
|    `-t`    | `--threads`       | `2`                                                      | number of threads to use for simulation                              |
|    `-v`    | `--vcf`           |                                                          | VCF-formatted file containing genotypes from which to create per-sample haplotypes |

### Output type
Mimick lets you specify different output fastq types regardless of the intended linked-read
simulation type. See [Data Formats](data_formats.md) for more information.

+++ Linked-read simulation
These are the options available specific to linked-read parameters, such as the average molecule length, etc.

{.compact}
| short name | long name             |                         default                          | description                                                                                                                                    | notes                                                                               |
|:----------:|:----------------------|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|
|            | `FASTA`               |                                                          | input fasta file(s)                                                                                                                            | REQUIRED                                                                            |
|    `-g`    | `--genomic-coverage`  |                           `30`                           | mean coverage target for simulated data                                                                                                        |                                                                                     |
|    `-i`    | `--insert-size`       |                          `500`                           | outer distance between the two ends in bp                                                                                                      | must be >`--read-lengths`                                                                 |
|    `-d`    | `--insert-stdev`      |                           `50`                           | standard deviation of --distance                                                                                                               |                                                                                     |
|    `-A`    | `--molecule-attempts` |                          `300`                           | number of attempts to create a molecule with <70% ambiguous bases before exiting with an error                                                 |                                                                                     |
|    `-C`    | `--molecule-coverage` |                          `0.2`                           | mean percent coverage per molecule if <1, else mean number of reads per molecule                                                               |                                                                                     |
|    `-L`    | `--molecule-length`   |                         `80000`                          | mean length of molecules in bp, drawn from exponential distribution                                                                            |                                                                                     |
|    `-N`    | `--molecules-per`     |                           `3`                            | mean number of unrelated molecules per barcode, drawn from an exponential distribution. If negative, (e.g. `-2`) will be fixed for that number |                                                                                     |
|    `-l`    | `--read-lengths`      |                        `150,150`                         | length of R1,R2 reads in bp                                                                                                                    | each must be >10 and separated by a comma, no spaces                                |
|    `-S`    | `--singletons`        |                           `0`                            | proportion of barcodes will only have a single read pair                                                                                       |

+++
