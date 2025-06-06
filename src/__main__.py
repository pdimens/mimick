#!/usr/bin/env python3

import glob
import os
import sys
import re
import math
import gzip
import multiprocessing
from itertools import product
import numpy as np
import rich_click as click
from rich.progress import Progress, TextColumn, TimeElapsedColumn, TaskProgressColumn, BarColumn
import pysam
from .classes import *
from .cli_classes import *
from .common import *
from .file_ops import *
from .simulate import *

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "Documentation: [link=https://pdimens.github.io/mimick/]https://pdimens.github.io/mimick/[/link]"
click.rich_click.OPTION_GROUPS = {
    "mimick": [
        {
            "name": "General Options",
            "options": ["--help", "--output-prefix", "--output-type", "--quiet", "--regions", "--seed", "--threads", "--version"],
            "panel_styles": {"border_style": "dim"}
        },
        {
            "name": "Read Simulation Parameters",
            "options": ["--coverage","--distance","--error","--extindels","--indels","--length","--mutation","--stdev"],            
            "panel_styles": {"border_style": "dim blue"}
        },
        {
            "name": "Linked Read Parameters",
            "options": ["--lr-type", "--molecule-coverage", "--molecule-length", "--molecules-per", "--singletons"],
            "panel_styles": {"border_style": "dim magenta"}
        },
    ]
}

@click.version_option("0.0.0", prog_name="mimick")
@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/SIM", show_default=True)
@click.option('-O','--output-type', help='output format of FASTQ files', type = click.Choice(["10x", "stlfr", "standard", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-q','--quiet', show_default = True, default = "0", type = click.Choice([0,1,2]), help = '`0` all output, `1` no progress bar, `2` no output')
@click.option('-r','--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('-t','--threads', help='number of threads to use for simulation', type=click.IntRange(min=1), default=2, show_default=True)
@click.option('-S','--seed', help='random seed for simulation', type=click.IntRange(min=1))
#Paired-end FASTQ simulation using pywgsim
@click.option('--coverage', help='mean coverage target for simulated data', show_default=True, default=30.0, type=click.FloatRange(min=0.05))
@click.option('--distance', help='outer distance between the two ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('--error', help='base error rate', default=0.02, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--extindels', help='indels extension rate', default=0.25, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--indels', help='indels rate', default=0.15, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--length', help='length of reads in bp', default=150, show_default=True, type=click.IntRange(min=30))
@click.option('--mutation', help='mutation rate', default=0.001, show_default=True, type=click.FloatRange(min=0))
@click.option('--stdev', help='standard deviation of --distance', default=50, show_default=True, type=click.IntRange(min=0))
#Linked-read simulation
@click.option('-l', '--lr-type', help='type of linked-read experiment', default = "haplotagging", show_default=True, show_choices=True, type= click.Choice(["10x", "stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-c','--molecule-coverage', help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-m','--molecule-length', help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=50))
@click.option('-n','--molecules-per', help='mean number of unrelated molecules per barcode per chromosome, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Normal distribution', default=2, show_default=True, type=int)
@click.option('-s','--singletons', help='proportion of barcodes will only have a single read pair', default=0, show_default=True, type=click.FloatRange(0,1))
@click.argument('barcodes', type = Barcodes())
@click.argument('fasta', type = click.Path(exists=True, dir_okay=False, resolve_path=True, readable=True), nargs = -1, required=True)
def mimick(barcodes, fasta, output_prefix, output_type, quiet, seed, regions, threads,coverage,distance,error,extindels,indels,length,mutation,stdev,lr_type, molecule_coverage, molecule_length, molecules_per, singletons):
    """
    Simulate linked-read FASTQ using genome haplotypes. Barcodes can be supplied one of two ways:
   
    1. let Mimick randomly generate barcodes based on a specification of `length,count`
        - two integers, comma-separated, no space
        - e.g. `16,400000` would generate 400,000 unique 16bp barcodes 
    2. you can provide a file of specific nucleotide barcodes, 1 per line

    You can specify the linked-read barcode chemistry to simulate via `--lr-type` as well as
    the output format of FASTQ files (default is the same as barcode type). For example, you
    can generate 96 barcodes (common haplotagging style), select `--lr-type stlfr`
    (combinatorial 3-barcode on R2 read), and have `--output-type tellseq` (`@seqid:barcode` header format).

    | --lr-type       | Format                                                 |
    |:----------------|:-------------------------------------------------------|
    | `10x`/`tellseq` | single barcode on R1                                   |
    | `haplotagging`  | R1 and R2 each have different combinatorial 2-barcodes |
    | `stlfr`         | combinatorial 3-barcode on R2                          |

    | --output-type  | Barcode Location                                      | Example                    |
    |:---------------|:------------------------------------------------------|:---------------------------|
    | `10x`          | start of R1 sequence                                  | `ATAGACCATAGA`GGACA...     |
    | `haplotagging` | sequence header as `BX:Z:ACBD`                        | `@SEQID BX:Z:A0C331B34D87` |
    | `standard`     | sequence header as `BX:Z:BARCODE`, no specific format | `@SEQID BX:Z:ATACGAGACA`   |
    | `stlfr`        | appended to sequence ID via `#1_2_3`                  | `@SEQID#1_354_39`          |
    | `tellseq`      | appended to sequence ID via `:ATCG`                   | `@SEQID:TATTAGCAC`         |
    """
    Container.CONSOLE = Console(stderr=True, log_path=False, quiet= quiet==2)
    Container.PROGRESS = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        transient=True,
        console=Container.CONSOLE,
        disable=quiet > 0
    )

    Container.OUT = os.path.dirname(output_prefix)
    os.makedirs(Container.OUT, exist_ok= True)
    Container.FASTADIR = fasta
    Container.PREFIX = os.path.basename(output_prefix)
    Container.RNG = np.random.default_rng(seed)
    Container.threads = threads
    Container.barcodetype=lr_type.lower()
    Container.coverage=coverage
    Container.error=error
    Container.distance=distance
    Container.stdev=stdev
    Container.length=length
    Container.mutation=mutation
    Container.indels=indels
    Container.extindels=extindels
    Container.mollen=molecule_length
    Container.singletons = singletons
    Container.used_bc = dict()
    if molecules_per == 0:
        error_terminate("The value for [yellow]--molecule-number[/] cannot be 0.")
    else:
        Container.molnum=molecules_per
    Container.molcov = molecule_coverage
    if isinstance(barcodes, str):
        Container.barcodepath=barcodes
    else:
        Container.barcodepath = f"{output_prefix}.generated.barcodes"
        bp,count = barcodes
        Container.CONSOLE.log(f'Generating {count} {bp}bp barcodes')
        with open(f"{Container.barcodepath}", "w") as bc_out:
            for i,bc in enumerate(generate_barcodes(bp),1):
                bc_out.write("".join(bc) + "\n")
                if i == count:
                    break
    if regions:
        # use the provided BED file
        Container.CONSOLE.log('Reading the input regions')
        intervals = readBED(regions)
    else:
        # derive intervals from the first FASTA file
        Container.CONSOLE.log(f'Inferring regions from [blue]{os.path.basename(fasta[0])}[/]')
        intervals = FASTAtoBED(fasta[0])
    if quiet == 0:
        Container.PROGRESS.start()
    pbar = Container.PROGRESS.add_task("[yellow]Progress", total=len(fasta) * len(intervals) + 1 + 1 + (2 * len(fasta)))
    Container.CONSOLE.log('Validating the input barcodes')
    try:
        with gzip.open(Container.barcodepath, 'rt') as filein:
            Container.barcodes, Container.barcodebp, Container.totalbarcodes = interpret_barcodes(filein, Container.barcodetype)
    except gzip.BadGzipFile:
        with open(Container.barcodepath, 'r') as filein:
            Container.barcodes, Container.barcodebp, Container.totalbarcodes = interpret_barcodes(filein, Container.barcodetype)
    except:
        error_terminate(f'Cannot open {Container.barcodepath} for reading')
    Container.PROGRESS.update(pbar, advance=1)
    # fill Container with wgsim and general linked-read parameters 
    Container.remainingbarcodes = Container.totalbarcodes
    Container.barcodeslist= open(f'{Container.OUT}/{Container.PREFIX}.barcodes', 'w')

    if Container.barcodetype in ["10x", "tellseq"]:
        # barcode at beginning of read 1
        Container.len_r1 = Container.length - Container.barcodebp
        if Container.len_r1 <= 5:
            error_terminate(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {Container.length}, Barcode length: {Container.barcodebp}')
        Container.len_r2 = Container.length
    elif Container.barcodetype == "stlfr":
        # barcode at the end of read 2
        Container.len_r1 = Container.length
        Container.len_r2 = Container.length - Container.barcodebp
        if Container.len_r2 <= 5:
            error_terminate(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {Container.length}, Barcode length: {Container.barcodebp}')
    else:
        # would be 4-segment haplotagging where AC on read 1 and BD on read 2, but in index read, so read not shortened
        Container.len_r1 = Container.length
        Container.len_r2 = Container.length

    if output_type:
        Container.outformat = output_type.lower()
    else:
        Container.outformat = Container.barcodetype

    if Container.outformat == "haplotagging":
        bc_range = [f"{i}".zfill(2) for i in range(1,97)]
        Container.bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
        # check that the haplotagging output format can support the supplied number of barcodes
        if Container.totalbarcodes > 96**4:
            error_terminate(f'The barcodes and barcode type supplied will generate a potenial {Container.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')

    if Container.outformat == "stlfr":
        bc_range = range(1, 1537)
        Container.bc_generator = product(bc_range, bc_range, bc_range)
        # check that the stlfr output format can support the supplied number of barcodes
        if Container.totalbarcodes > 1537**3:
            error_terminate(f'The barcodes and barcode type supplied will generate a potenial {Container.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')

    Container.ffiles = fasta 
    Container.regioncoverage = Container.coverage/len(Container.ffiles)

    for k,s in enumerate(Container.ffiles,1):
        Container.CONSOLE.rule(f"[bold]Haplotype {k}: {os.path.basename(s)}", style = "blue")
        Container.hapnumber = f'{k}'
        try:
            Container.CONSOLE.log(f"Indexing FASTA file")
            Container.ffile = os.path.abspath(s)
            pysam.faidx(Container.ffile)
            if Container.ffile.lower().endswith(".gz") and not os.path.exists(Container.ffile + ".gzi"):
                os.system(f"bgzip --reindex {Container.ffile}")
        except Exception as e:
            error_terminate(f'Failed to index {Container.ffile}. Error reported by samtools:\n{e}', False)
        for w in intervals:
            regiontext = f"{w.chrom}:{w.start}-{w.end}"
            rule = "─" * int((53 - len(regiontext)) / 2)
            Container.CONSOLE.log(f"[green]{rule} [default]{regiontext} [green]{rule}[/]")
            try:
                LinkedSim(w,Container)
            except KeyboardInterrupt:
                Container.PROGRESS.stop()
                mimick_keyboardterminate()
            Container.PROGRESS.update(pbar, advance=1)

    Container.barcodeslist.close()
    # gzip multiprocessing
    allfastq = glob.glob(os.path.abspath(Container.OUT) + '/*.fq')
    chunk_size = len(allfastq)/Container.threads
    slices = Chunks(allfastq, math.ceil(chunk_size))
    processes = []
    Container.CONSOLE.rule("Finalizing Outputs", style = "blue")
    for i,sli in enumerate(slices):
        p = multiprocessing.Process(target=BGzipper, args=(sli,))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
        Container.PROGRESS.update(pbar, advance=1)
    tb = log_table()
    tb.add_row("Concatenating and gzipping", "WGSim GFF files")
    Container.CONSOLE.log(tb)
    with gzip.open(f'{Container.OUT}/{Container.PREFIX}.mutations.gff.gz', "wb", compresslevel=6) as gff:
        gff.write("##gff-version 3\n#\n".encode("utf-8"))
        for i in range(1, Container.threads + 1):
            with open(f'{Container.OUT}/{Container.PREFIX}.p{i}.wgsim.mutations', 'r') as mut:
                for line in mut:
                    if not line.startswith("#"):
                        gff.write(line.encode("utf-8"))
            os.remove(f'{Container.OUT}/{Container.PREFIX}.p{i}.wgsim.mutations')
    Container.PROGRESS.update(pbar, completed=100)
    Container.PROGRESS.stop()
    Container.CONSOLE.rule("Finished Successfully", style = "green bold")

if __name__ =='__main__':
    try:
        mimick()
    except KeyboardInterrupt:
        mimick_keyboardterminate()

