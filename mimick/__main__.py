#!/usr/bin/env python3

import os
import sys
import rich_click as click
from .cli_classes import ReadLengths

config = click.RichHelpConfiguration(
    max_width=80,
    theme = "magenta2-slim",
    use_markdown=True,
    show_arguments=False,
    style_options_panel_border = "magenta",
    style_commands_panel_border = "yellow",
    style_option_default= "dim",
    style_deprecated="dim red",
    options_table_column_types = ["opt_long", "opt_short", "help"],
    options_table_help_sections = ["required", "help", "default"]
)

@click.version_option("0.0.0", prog_name="mimick")
@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
@click.rich_config(config)
@click.option('-c', '--circular', panel = "General Options", is_flag= True, default = False, help = 'contigs are circular/prokaryotic')
@click.option('-o', '--output-prefix', panel = "General Options", help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/", show_default=True)
@click.option('-q', '--quiet', panel = "General Options", is_flag= True, default = False, help = 'toggle to hide progress bar')
@click.option('-t', '--threads', panel = "General Options", help='number of threads to use for multi-sample simulation', type=click.IntRange(min=1), default=2, show_default=True)
@click.option('-s', '--seed', panel = "General Options", help='random seed for simulation', type=click.IntRange(min=-1, clamp = True), default=-1)
@click.option('-f', '--format', 'fmt', panel = "Linked-Read Simulation", help='FASTQ output format',show_default=True, default = "standard:haplotagging", type = click.Choice(["10x", "stlfr", "standard", "standard:haplotagging", "standard:stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-g', '--genomic-coverage', panel = "Linked-Read Simulation", help='mean coverage (depth) target for simulated data', show_default=True, default=30.0, type=click.FloatRange(min=0.05))
@click.option('-i', '--insert-size', panel = "Linked-Read Simulation", help='outer distance between the two read ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('-d', '--insert-stdev', panel = "Linked-Read Simulation", help='standard deviation for `--insert-size`', default=50, show_default=True, type=click.IntRange(min=0))
@click.option('-A', '--molecule-attempts', panel = "Linked-Read Simulation", help='how many tries to create a molecule with <70% ambiguous bases', show_default=True, default=300, type=click.IntRange(min=5))
@click.option('-C', '--molecule-coverage', panel = "Linked-Read Simulation", help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-L', '--molecule-length', panel = "Linked-Read Simulation", help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=650))
@click.option('-N', '--molecules-per', panel = "Linked-Read Simulation", help='mean number of unrelated molecules per barcode per chromosome, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a distribution', default=2, show_default=True, type=int)
@click.option('-l', '--read-lengths', panel = "Linked-Read Simulation", help='length of R1,R2 sequences in bp', default="150,150", show_default=True, type=ReadLengths())
@click.option('-S', '--singletons', panel = "Linked-Read Simulation", help='proportion of barcodes that will only have one read pair', default=0, show_default=True, type=click.FloatRange(0,1))
@click.option('-v', '--vcf', help='VCF-formatted file containing genotypes from which to create per-sample haplotypes', type = click.Path(exists=True, dir_okay=False, resolve_path=True, readable=True))
@click.argument('fasta', type = click.Path(exists=True, dir_okay=False, resolve_path=True, readable=True), nargs = -1, required=True)
def mimick(fasta, circular, quiet, output_prefix, fmt, seed, threads,genomic_coverage,insert_size,read_lengths,insert_stdev, molecule_coverage, molecule_attempts, molecule_length, molecules_per, singletons, vcf):
    """
    Simulate linked-read FASTQ data for one or many individuals

    There are two modes of operation:
    1. Input one or more **FASTA** files (haplotypes) to simulate linked reads for a single individual.
    2. Input one **FASTA** and **VCF** file to simulate linked reads for all samples in the VCF file with haplotypes reflective of their SNPs and indels.

    With the exception of `10x`, all other formats are demultiplexed. Output can be in `standard:chemistry` format (e.g. `standard:stlfr` outputs standard
    format with stLFR-style barcodes), where the barcode is encoded as a `BX:Z:` tag and a `VX:i` validation tag. Below are the common linked-read chemistries
    (to be used in `--format`) and their configurations:

    | chemistry    | `--read-lengths` | Description          | FASTQ format          |
    |:-------------|:----------------:|:---------------------|:----------------------|
    | 10x          |    `134,150`     | single barcode on R1 | barcode inline on R1  |
    | haplotagging |    `150,150`     | 2-barcodes in I1/I2  | BX:Z:AxxCxxBxxDxx tag |
    | stlfr        |    `150,108`     | 3-barcode on R2      | @SEQID#1_2_3          |
    | tellseq      |    `132,150`     | single barcode on R1 | @SEQID:ATGC           |
    """
    if fmt == "10x":
        fmt = "tenx"
    os.environ['PYTHON_JULIACALL_THREADS'] = f"{threads}"
    os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS'] = "yes"
    from juliacall import Main as jl
    jl.seval("using MimickLinkedReads")
    if vcf:
        jl.mimick(
            fasta[0],
            vcf,
            fmt,
            outdir = output_prefix,
            coverage = genomic_coverage,
            n_molecules = molecules_per,
            mol_coverage = molecule_coverage,
            mol_length = molecule_length,
            insert_length = insert_size,
            insert_stdev = insert_stdev,
            read_length = jl.collect(jl.Int, read_lengths),
            singletons = singletons,
            circular = circular,
            attempts = molecule_attempts,
            seed = seed,
            quiet = quiet
    )
    else:
        jl.mimick(
            jl.collect(fasta),
            fmt,
            prefix = output_prefix,
            coverage = genomic_coverage,
            n_molecules = molecules_per,
            mol_coverage = molecule_coverage,
            mol_length = molecule_length,
            insert_length = insert_size,
            insert_stdev = insert_stdev,
            read_length = jl.collect(jl.Int, read_lengths),
            singletons = singletons,
            circular = circular,
            attempts = molecule_attempts,
            seed = seed,
            quiet = quiet
        )

def mimick_test():
    """
    A simple function to make sure JuliaCall is working properly after a conda installation
    """
    try:
        from juliacall import Main as jl
        jl.seval("using MimickLinkedReads")
        print("Success!")
    except Exception as e:
        print(e)
        print("Failure!")

if __name__ =='__main__':
    try:
        mimick()
    except KeyboardInterrupt:
        exit(1)
    except Exception as e:
        print(e)
        sys.exit(1)
