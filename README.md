# Mimick
Linked-read sequence simulator

![mimick_logo](/docs/static/mimick_logo.png)

[![gitHub release](https://img.shields.io/github/v/release/pdimens/mimick?style=for-the-badge&logo=anaconda&logoColor=ffffff)](https://github.com/pdimens/mimick/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-docs-daa355?style=for-the-badge&logo=quicklook&logoColor=ffffff)](https://pdimens.github.io/mimick)
[![harpy badge](https://custom-icon-badges.demolab.com/badge/-Harpy-79a9b9?style=for-the-badge&logo=package&logoColor=ffffff)](https://www.github.com/pdimens/harpy)
[![visor badge](https://custom-icon-badges.demolab.com/badge/-VISOR-12922e?style=for-the-badge&logo=package&logoColor=ffffff)](https://github.com/davidebolo1993/VISOR)

Originally known as XENIA from the [VISOR](https://github.com/davidebolo1993/VISOR) project, Mimick is a 
simulator for linked-read FASTQ data. Mimick allows you to simulate an
arbitrary number of haplotypes, set overall coverage, molecule coverage,
and choose what kind of linked reads you want.

## Supported Linked-Read Types
- 10X
- Haplotagging
- stLFR
- TELLseq

## Simulation parameters
- output FASTQ format
- overall coverage depth
- average molecule length
- molecule coverage / reads per molecule
- molecules per barcode (barcode convolution)
- proportion of singletons (unlinked barcodes)
- standard Illumina read characteristics e.g. read length, insert size, etc.

## Standout Features
Other than the fun name and logo, Mimick is an improvement over existing linked-read simulators in multiple ways:

1. It's the only simulator (we are aware of) that isn't configured for discontinued-in-2019 10X linked-read chemistry and is instead 
generalized for existing options, both in terms of data formats and the simulation process itself.
2. Circular DNA support. Yay prokaryotes!
3. Mimick provides more parameters to tune your simulations for realistic linked-read library simulation in the form of singletons and 
molecule coverage. These characteristics are **very** important regarding the performance of a linked-read library.
4. As of version 2.0, Mimick uses a barcode-first simulation approach, which allows barcodes to be shared **across**
chromosomes/contigs **and** haplotypes. This form of barcode sharing is a common phenomenon in real linked-read
libraries, but a characteristic existing simulators don't capture (e.g. XENIA only allowed barcode sharing within
a chromosome within a haplotype). The documentation explains this in better detail.
5. As of version 3.0 (upcoming), it supports multi-sample simulation by way of one FASTA and one VCF as input. Sample haplotypes
are made by applying SNP and indel variants from VCF to the contigs in the FASTA.
6. It's fast. The Julia version (v3+) is a signficant speedup for single-sample simulation and the multi-sample simulation
is parallelized across samples.

### Authors

<img src="https://avatars.githubusercontent.com/u/19176506?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@pdimens](https://github.com/pdimens) (Mimick)

<img src="https://avatars.githubusercontent.com/u/39052119?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@davidebolo1993](https://github.com/davidebolo1993) (VISOR)

> [!NOTE]
> Why name it "mimick"? Well, this software mimics linked-read data, I have an affinity for naming software after
> [fictional monsters](https://en.wikipedia.org/wiki/Mimic_(Dungeons_%26_Dragons)) and "mimick" (with a "k") is the old-English
> spelling of the word, leaving `mimic` available for some other bioinformatician to use for a less farcical reason. Despite the
> lore of mimics being deadly traps, this software is anything but, I promise.

