# Using Mimick

```using mimick
mimick options... BARCODES FASTA1 FAST2...
```
Use `--help` or `-h` or call `mimick` without arguments to call up the docstring.

The minimum input files required by Mimick is a single FASTA file, uncompressed or **bgzip compressed**.
The result will be a single set of paired-end reads, a GFF file of mutations, and a manifest of all the molecules created.
Mimick scales with the number of threads provided.

If you want to separate out the haplotypes from the final output, you can
leverage the fact that all the simulated reads start with `@HAP:X_` where `X` is the haplotype number, starting with `1`,
corresponding to the order in which FASTA files were provided. The read names also
include source contig names and other identifying features that can be extracted similarly.

```bash
for i in {1..2}; do
    zgrep -A3 \"^@HAP:X_\" output_prefix.R$i.fq.gz | gzip > out.R$i.fq.gz
done
```

<!-- tabs:start -->

#### **Required arguments**
Running Mimick will always require a barcode specification and at least one input FASTA file.

### Barcodes
#### Randomly generate
Mimick lets you put in length and count parameters, which it will use to randomly generate barcodes.
The format is `length,count` (no spaces), where `length` is the base-pair length the barcodes should be and
`count` is how many barcodes it should generate of length `length`. For example, if you specify `16,4000000`,
Mimick will generate 4 million unique 16bp barcodes, effectively mimicking 10X barcodes (see what I did there?).
Mimick will write a file containing the barcodes it generated. In practice, this would look something like:
```randomly generated barcodes
mimick --lr-type 10x 16,4000000 hap1.fasta hap2.fasta
#                 16bp^ ^4 million barcodes
```

#### Specific barcodes
Alternatively, if you have a set of barcodes you absolutely want to use, just put the filename as the first positional argument.
In practice, this would look something like:
```barcodes as a file
mimick --lr-type 10x barcodes.txt hap1.fasta hap2.fasta
#                    ^file of barcodes
```

### FASTA file(s)
You will need at least 1 fasta file as input, which goes at the very end of the command. If providing more than 1
FASTA for a non-haploid species, it's assumed that each fasta file is a different haplotype **of the same genome**.
There is no strict enforcement of the fasta files being from the same genome, it's just how you'll probably want to use
the simulator (but I'm open to being surprised). 

```fasta inputs
mimick --lr-type 10x 16,4000000 hap1.fasta hap2.fasta
#                    ^barcodes  ^fasta 1   ^fasta 2
```
#### Circular DNA
The `--circular` toggle/flag tells Mimick to treat each contig
within each FASTA as circular when creating molecules-- this is probably how you'll want to simulate prokaryotic/microbial
genomes. When using `--circular`, you'll see in the `.molecules` file that some molecule end positions may occur before the start
positions, which is an artifact of the circularization. Those positions reflect the start/end positions on the linear sequence in
the FASTA file and should be interpreted as "started at `start` and reached the end of the contig, then wrapped around to the
beginning of the contig and kept going until `end`." The corresponding length of the molecule will be accurate and the math
should make sense as well: $end = start + molecule\_size - contig\_size$

#### **Detailed Options**

### General options
These options control inputs/outputs and resources

| short name | long name         | default         | description                                                          |
|:----------:|:------------------|:----------------|:---------------------------------------------------------------------|
|    `-C`    | `--circular`      | `False`         | toggle to let Mimick know the input contigs are circular/prokaryotic |
|    `-o`    | `--output-prefix` | `simulated/SIM` | output file prefix                                                   |
|    `-O`    | `--output-type`   | same as input   | output format of FASTQ files                                         |
|    `-q`    | `--quiet`         | `0`             | `0` all output, `1` no progress bar, `2` no output                   |
|    `-r`    | `--regions`       |                 | one or more regions to simulate, in BED format                       |
|    `-S`    | `--seed`          |                 | random seed for simulation                                           |
|    `-t`    | `--threads`       | `2`             | number of threads to use for simulation                              |

#### Output type
Mimick lets you specify different output fastq types regardless of the intended linked-read
simulation type. See [Data Formats](data_formats.md) for more information.

### FASTQ simulation with pywgsim
These options govern how `wgsim` will simulate FASTQ files from genomic regions. These are no short names for
these options.

|     name      | default | description                               | notes                                             |
|:-------------:|:--------|:------------------------------------------|:--------------------------------------------------|
| `--coverage`  | `30`    | mean coverage target for simulated data   |                                                   |
| `--distance`  | `500`   | outer distance between the two ends in bp | must be >`--length`                               |
|   `--error`   | `0.02`  | base error rate                           | must be between 0-1, will be fixed for this value |
| `--extindels` | `0.25`  | indels extension rate                     | must be between 0-1                               |
|  `--indels`   | `0.15`  | indels creation rate                      | must be between 0-1                               |
|  `--length`   | `150`   | length of reads in bp                     | must be >30                                       |
| `--mutation`  | `0.001` | mutation rate                             | must be between 0-1                               |
|   `--stdev`   | `50`    | standard deviation of --distance          |                                                   |

### Linked-read simulation
| short name | long name             | default        | description                                                                                                                                                                                |
|:----------:|:----------------------|:---------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|    `-l`    | `--lr-type`           | `haplotagging` | type of linked-read experiment                                                                                                                                                             |
|    `-a`    | `--molecule-attempts` | `300`          | number of attempts to create a molecule with <70% ambiguous bases before exiting with an error                                                                                             |
|    `-c`    | `--molecule-coverage` | `0.2`          | mean percent coverage per molecule if <1, else mean number of reads per molecule'                                                                                                          |
|    `-m`    | `--molecule-length`   | `80000`        | mean length of molecules in bp'                                                                                                                                                            |
|    `-n`    | `--molecules-per`     | `3`            | mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Normal distribution |
|    `-s`    | `--singletons`        | `0`            | proportion of barcodes will only have a single read pair                                                                                                                                   |

#### **All Options**
For completeness, the table below is all the command line arguments and options

| short name | long name             | default         | description                                                                                                                                    | notes                                                                                                   |
|:----------:|:----------------------|:----------------|:-----------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|
|            | `BARCODES`            |                 | input barcode file or length,count                                                                                                             | REQUIRED                                                                                                |
|            | `FASTA`               |                 | input fasta file(s)                                                                                                                            | REQUIRED                                                                                                |
|    `-C`    | `--circular`          | `False`         | toggle to let Mimick know the input contigs are circular/prokaryotic                                                                           | cannot be used with `-r/--regions`                                                                      |
|    `-o`    | `--output-prefix`     | `simulated/SIM` | output file prefix                                                                                                                             |                                                                                                         |
|    `-O`    | `--output-type`       | `standard`      | output format of FASTQ files                                                                                                                   |                                                                                                         |
|    `-q`    | `--quiet`             | `0`             | `0` all output, `1` no progress bar, `2` no output                                                                                             |                                                                                                         |
|    `-r`    | `--regions`           |                 | one or more regions to simulate, in BED format                                                                                                 |                                                                                                         |
|    `-S`    | `--seed`              |                 | random seed for simulation                                                                                                                     | optional and useful for reproducibility                                                                 |
|    `-t`    | `--threads`           | `2`             | number of threads to use for simulation                                                                                                        |                                                                                                         |
|            | `--coverage`          | `30`            | mean coverage target for simulated data                                                                                                        |                                                                                                         |
|            | `--distance`          | `500`           | outer distance between the two ends in bp                                                                                                      | must be >`--length`                                                                                     |
|            | `--error`             | `0.02`          | base error rate                                                                                                                                | must be between 0-1, will be fixed for this value                                                       |
|            | `--extindels`         | `0.25`          | indels extension rate                                                                                                                          | must be between 0-1                                                                                     |
|            | `--indels`            | `0.15`          | indels creation rate                                                                                                                           | must be between 0-1                                                                                     |
|            | `--length`            | `150`           | length of reads in bp                                                                                                                          | must be >30                                                                                             |
|            | `--mutation`          | `0.001`         | mutation rate                                                                                                                                  | must be between 0-1                                                                                     |
|            | `--stdev`             | `50`            | standard deviation of --distance                                                                                                               |                                                                                                         |
|    `-l`    | `--lr-type`           | `haplotagging`  | type of linked-read experiment                                                                                                                 | determines if barcodes are combinatorial, how, and where read lengths get shortened to accommodate them |
|    `-a`    | `--molecule-attempts` | `300`           | number of attempts to create a molecule with <70% ambiguous bases before exiting with an error                                                 |                                                                                                         |
|    `-c`    | `--molecule-coverage` | `0.2`           | mean percent coverage per molecule if <1, else mean number of reads per molecule                                                               |                                                                                                         |
|    `-m`    | `--molecule-length`   | `80000`         | mean length of molecules in bp, drawn from exponential distribution                                                                            |                                                                                                         |
|    `-n`    | `--molecules-per`     | `3`             | mean number of unrelated molecules per barcode, drawn from an exponential distribution. If negative, (e.g. `-2`) will be fixed for that number |                                                                                                         |
|    `-s`    | `--singletons`        | `0`             | proportion of barcodes will only have a single read pair                                                                                       |

<!-- tabs:end -->
