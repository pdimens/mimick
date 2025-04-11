# Using Mimick

```using mimick
mimick options... BARCODES FASTA1 FAST2...
```
Use `--help` or `-h` or call `mimick` without arguments to call up the docstring.


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
You will need at least 1 fasta file as input, which goes at the very end of the command. It's assumed
that each fasta file is a different haplotype **of the same genome**.
```fasta inputs
mimick --lr-type 10x 16,4000000 hap1.fasta hap2.fasta
#                    ^barcodes  ^fasta 1   ^fasta 2
```

#### **Detailed Options**

### General options
These options control inputs/outputs and resources

| short name | long name | description | default |
|:---:|:------|:-------------|:--------|
`-o` |`--output-prefix` | output file prefix | `simulated/SIM`|
`-O` |`--output-format` | output format of FASTQ files | same as input |
`-r` |`--regions` | one or more regions to simulate, in BED format | |
`-t` |`--threads` | number of threads to use for simulation | `2` |

#### Output format
Mimick lets you specify different output fastq types regardless of the intended linked-read
simulation type. See [Data Formats](data_formats.md) for more information.

### FASTQ simulation with pywgsim
These options govern how `wgsim` will simulate FASTQ files from genomic regions. These are no short names for
these options.

| name | description | default | notes |
|:---:|:------|:-------------|:--------|
| `--coverage` | 'mean coverage target for simulated data' | `30` | |
| `--distance` | 'outer distance between the two ends in bp' | `500` | must be >`--length` |
| `--error` | 'base error rate' | `0.02`| must be between 0-1, will be fixed for this value |
| `--extindels` | 'indels extension rate' | `0.25` | must be between 0-1 |
| `--indels` | 'indels creation rate' | `0.15` | must be between 0-1 |
| `--length` | 'length of reads in bp' | `150` | must be >30 |
| `--mutation` | 'mutation rate' | `0.001` | must be between 0-1 |
| `--stdev` | 'standard deviation of --distance' | `50` | |

### Linked-read simulation
| short name | long name | description | default |
|:---:|:------|:-------------|:--------|
|`-l` | `--lr-type` | type of linked-read experiment | `haplotagging` |
|`-c` | `--molecule-coverage` | mean percent coverage per molecule if <1, else mean number of reads per molecule' | `0.2` |
|`-m` | `--molecule-length` | mean length of molecules in bp' | `80000` |
|`-n` | `--molecule-number` | mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Poisson distribution | 3 |

#### **All Options**
For completeness, the table below is all the command line arguments and options

| short name | long name | description | default | notes |
|:---:|:------|:-------------|:--------|:----------------|
|     | `BARCODES` | input barcode file or length,count | | REQUIRED |
|     | `FASTA` | input fasta file(s) | | REQUIRED |
|`-o` |`--output-prefix` | output file prefix | `simulated/SIM`| |
|`-O` |`--output-format` | output format of FASTQ files | `standard` | |
|`-r` |`--regions` | one or more regions to simulate, in BED format | | |
|`-t` |`--threads` | number of threads to use for simulation | `2` | |
|     | `--coverage` | 'mean coverage target for simulated data' | `30` | |
|     | `--distance` | 'outer distance between the two ends in bp' | `500` | must be >`--length` |
|     | `--error` | 'base error rate' | `0.02`| must be between 0-1, will be fixed for this value |
|     | `--extindels` | 'indels extension rate' | `0.25` | must be between 0-1 |
|     | `--indels` | 'indels creation rate' | `0.15` | must be between 0-1 |
|     | `--length` | 'length of reads in bp' | `150` | must be >30 |
|     | `--mutation` | 'mutation rate' | `0.001` | must be between 0-1 |
|     | `--stdev` | 'standard deviation of --distance' | `50` | |
|`-l` | `--lr-type` | type of linked-read experiment | `haplotagging` | |
|`-c` | `--molecule-coverage` | mean percent coverage per molecule if <1, else mean number of reads per molecule' | `0.2` | |
|`-m` | `--molecule-length` | mean length of molecules in bp' | `80000` | |
|`-n` | `--molecule-number` | mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Poisson distribution | `3`| |

<!-- tabs:end -->
