#!/usr/bin/env python3

import shutil
import os
import rich_click as click
from rich import print as rprint
from .cli_classes import ReadLengths

os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS'] = "yes"

if shutil.which("julia"):
    os.environ['PYTHON_JULIAPKG_EXE'] = shutil.which("julia")

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
@click.option_panel(
    "General Options",
    options = ["--circular", "--output-prefix", "--quiet", "threads", "--seed", "--vcf", "--version", "--help"],
    title_style="yellow" 
)
@click.option_panel(
    "Linked-Read Simulation",
    options = ["--format", "--genomic-coverage", "--insert-size", "--insert-stdev", "--molecule-attempts", "--molecule-coverage", "--molecule-length", "--molecules-per", "--read-lengths", "--singletons"],
    title_style="yellow"
)
@click.option('-c', '--circular', is_flag= True, default = False, help = 'contigs are circular/prokaryotic')
@click.option('-o', '--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/", show_default=True)
@click.option('-q', '--quiet', is_flag= True, default = False, help = 'toggle to hide progress bar')
@click.option('-t', '--threads', help='number of threads to use for multi-sample simulation', type=click.IntRange(min=2), default=2, show_default=True)
@click.option('-s', '--seed', help='random seed for simulation', type=click.IntRange(min=-1, clamp = True), default=-1)
@click.option('-f', '--format', 'fmt', help='FASTQ output format',show_default=True, default = "standard:haplotagging", type = click.Choice(["10x", "stlfr", "standard", "standard:haplotagging", "standard:stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-g', '--genomic-coverage', help='mean coverage (depth) target for simulated data', show_default=True, default=30.0, type=click.FloatRange(min=0.05))
@click.option('-i', '--insert-size', help='outer distance between the two read ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('-d', '--insert-stdev', help='standard deviation for `--insert-size`', default=50, show_default=True, type=click.IntRange(min=0))
@click.option('-A', '--molecule-attempts', help='how many tries to create a molecule with <70% ambiguous bases', show_default=True, default=300, type=click.IntRange(min=5))
@click.option('-C', '--molecule-coverage', help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-L', '--molecule-length', help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=650))
@click.option('-N', '--molecules-per', help='mean number of unrelated molecules per barcode per chromosome, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a distribution', default=2, show_default=True, type=int)
@click.option('-l', '--read-lengths', help='length of R1,R2 sequences in bp', default="150,150", show_default=True, type=ReadLengths())
@click.option('-S', '--singletons', help='proportion of barcodes that will only have one read pair', default=0, show_default=True, type=click.FloatRange(0,1))
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
    os.environ['PYTHON_JULIACALL_THREADS'] = f"{threads}"
    from juliacall import Main as jl
    try:
        jl.seval("using MimickLinkedReads")
    except:
        mimick_install()

    fmt = "tenx" if fmt == "10x" else fmt

    try:
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
    except KeyboardInterrupt:
        rprint("[yellow]Terminating Mimick")
    except Exception as e:
        rprint(f"[red]{e}")

def mimick_install():
    """Install the Mimick Julia backend"""
    os.environ['PYTHON_JULIACALL_THREADS'] = "2"
    from juliacall import Main as jl
    import mimick
    rprint("[magenta]Installing the MimickLinkedreads Julia backend. This typically only needs to happen once.")
    try:
        jl.seval('using Pkg')
        #jl.Pkg.develop(path=os.path.join(os.path.dirname(mimick.__path__[0]), "MimickLinkedReads.jl/"))
        jl.Pkg.develop(path=os.path.join(mimick.__path__[0], "MimickLinkedReads.jl/"))
        jl.seval('using MimickLinkedReads')
    except Exception as e:
        print("Failed to install MimickLinkedReads.jl. See the Julia error:")
        print(e)

def mimick_test():
    """
    A simple function to make sure JuliaCall is working properly after a conda installation
    """
    from juliacall import Main as jl
    try:
        jl.seval("using MimickLinkedReads")
        print("Success!")
    except Exception as e:
        print(e)
        print("Failure!")

