#!/usr/bin/env python3

import glob
import os
import sys
import rich_click as click
import re
import math
import gzip
import multiprocessing
from itertools import product
from .classes import *
from .cli_classes import *
from .common import *
from .file_ops import *
from .simulate import *
import pysam

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
            "options": ["--help", "--output-prefix", "--output-type", "--regions", "--threads", "--version"],
            "panel_styles": {"border_style": "dim"}
        },
        {
            "name": "Read Simulation Parameters",
            "options": ["--coverage","--distance","--error","--extindels","--indels","--length","--mutation","--stdev"],            
            "panel_styles": {"border_style": "dim blue"}
        },
        {
            "name": "Linked Read Parameters",
            "options": ["--lr-type", "--molecule-coverage", "--molecule-length", "--molecule-number"],
            "panel_styles": {"border_style": "dim magenta"}
        },
    ]
}

@click.version_option("0.0.0", prog_name="mimick")
@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/SIM", show_default=True)
@click.option('-O','--output-type', help='output format of FASTQ files', type = click.Choice(["10x", "stlfr", "standard", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-r','--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('-t','--threads', help='number of threads to use for simulation', type=click.IntRange(min=1), default=2, show_default=True)
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
@click.option('-n','--molecule-number', help='mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Poisson distribution', default=3, show_default=True, type=int)
@click.argument('barcodes', type = Barcodes())
@click.argument('fasta', type = click.Path(exists=True, dir_okay=False, resolve_path=True, readable=True), nargs = -1, required=True)
def mimick(barcodes, fasta, output_prefix, output_type, regions, threads,coverage,distance,error,extindels,indels,length,mutation,stdev,lr_type, molecule_coverage, molecule_length, molecule_number):
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

    | --lr-type | Format |
    |:------------------|:-------|
    |`10x`/`tellseq`   | single barcode on R1 |
    |`haplotagging`  | R1 and R2 each have different combinatorial 2-barcodes |
    |`stlfr`         | combinatorial 3-barcode on R2 |

    | --output-type | Barcode Location | Example |
    |:-----------------|:-------|:---------------------|
    |`10x`           | start of R1 sequence | `ATAGACCATAGA`GGACA... |
    |`haplotagging`  | sequence header as `BX:Z:ACBD` |  `@SEQID BX:Z:A0C331B34D87` |
    |`standard`      | sequence header as `BX:Z:BARCODE`, no specific format | `@SEQID BX:Z:ATACGAGACA` |
    |`stlfr`         | appended to sequence ID via `#1_2_3` | `@SEQID#1_354_39` |
    |`tellseq`       | appended to sequence ID via `:ATCG` | `@SEQID:TATTAGCAC` |
    """
    c.OUT = os.path.dirname(output_prefix)
    os.makedirs(c.OUT, exist_ok= True)
    c.FASTADIR = fasta
    c.PREFIX = os.path.basename(output_prefix)
    c.threads = threads
    c.barcodetype=lr_type.lower()
    c.coverage=coverage
    c.error=error
    c.distance=distance
    c.stdev=stdev
    c.length=length
    c.mutation=mutation
    c.indels=indels
    c.extindels=extindels
    c.molnum=molecule_number
    c.mollen=molecule_length
    c.molcov=molecule_coverage
    if isinstance(barcodes, str):
        c.barcodepath=barcodes
    else:
        c.barcodepath = f"{output_prefix}.generated.barcodes"
        bp,count = barcodes
        mimick_console.log(f'Generating {count} {bp}bp barcodes')
        with open(f"{c.barcodepath}", "w") as bc_out:
            for i,bc in enumerate(generate_barcodes(bp),1):
                bc_out.write("".join(bc) + "\n")
                if i == count:
                    break
    if regions:
        # use the provided BED file
        mimick_console.log('Reading the input regions')
        intervals = readBED(regions)
    else:
        # derive intervals from the first FASTA file
        mimick_console.log('Inferring regions from the first FASTA file')
        intervals = FASTAtoBED(fasta[0])
    
    mimick_console.log('Validating the input barcodes')
    try:
        with gzip.open(c.barcodepath, 'rt') as filein:
            c.barcodes, c.barcodebp, c.totalbarcodes = interpret_barcodes(filein, c.barcodetype)
    except gzip.BadGzipFile:
        with open(c.barcodepath, 'r') as filein:
            c.barcodes, c.barcodebp, c.totalbarcodes = interpret_barcodes(filein, c.barcodetype)
    except:
        click.BadParameter(f'Cannot open {c.barcodepath} for reading', param_hint= "BARCODES")
        #mimick_errorterminate(f'[Error] Cannot open {c.barcodepath} for reading', False)
    
    # fill c with wgsim and general linked-read parameters 
    c.remainingbarcodes = c.totalbarcodes
    c.barcodelist= open(f'{c.OUT}/{c.PREFIX}.barcodes', 'w')

    if c.barcodetype in ["10x", "tellseq"]:
        # barcode at beginning of read 1
        c.len_r1 = c.length - c.barcodebp
        if c.len_r1 <= 5:
            click.BadParameter(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', param_hint='--lr-type BARCODES')
            #mimick_errorterminate(f'[Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', False)
        c.len_r2 = c.length
    elif c.barcodetype == "stlfr":
        # barcode at the end of read 2
        c.len_r1 = c.length
        c.len_r2 = c.length - c.barcodebp
        if c.len_r2 <= 5:
            click.BadParameter(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', param_hint='--lr-type BARCODES')
            #mimick_errorterminate(f'[Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', False)
    else:
        # would be 4-segment haplotagging where AC on read 1 and BD on read 2
        c.len_r1 = c.length - c.barcodebp
        c.len_r2 = c.length - c.barcodebp
        if c.len_r1 <= 5 or c.len_r2 <= 5:
            click.BadParameter(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', param_hint='--lr-type BARCODES')
            #mimick_console.log(f'[Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', False)
            sys.exit(1)
    if output_type:
        c.outformat = output_type.lower()
    else:
        c.outformat = c.barcodetype
    if c.outformat == "haplotagging":
        bc_range = [f"{i}".zfill(2) for i in range(1,97)]
        c.bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
    if c.outformat == "stlfr":
        bc_range = range(1, 1537)
        c.bc_generator = product(bc_range, bc_range, bc_range)
    c.ffiles = fasta 
    c.regioncoverage = c.coverage/len(c.ffiles)
    
    # check that the haplotagging output format can support the supplied number of barcodes
    if c.outformat == "haplotagging":
        if c.totalbarcodes > 96**4:
            click.BadParameter('The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes', param_hint= '--lr-type --output-format')
            #mimick_errorterminate(f'[Error] The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')
    # check that the stlfr output format can support the supplied number of barcodes
    if c.outformat == "stlfr":
        if c.totalbarcodes > 1537**3:
            click.BadParameter('The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {1537**3} barcodes', param_hint= '--lr-type --output-format')
            #mimick_errorterminate(f'[Error] The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')

    for k,s in enumerate(c.ffiles):
        haplotype_table = log_table()
        haplotype_table.add_row(f'Processing haplotype {k+1}', f'[blue]{os.path.basename(s)}[/]')
        mimick_console.log(haplotype_table)
        #mimick_console.log(f'Processing haplotype {k+1}: [blue]{os.path.basename(s)}[/]')
        c.hapnumber = f'{k+1}'
        c.ffile = c.ffiles[k]
        for w in intervals:
            region_table = log_table()
            region_table.add_row('Region being simulated', f'[green]{w.chrom}:{w.start}-{w.end}[/]')
            mimick_console.log(region_table)
            try:
                LinkedSim(w,c)
                
            except KeyboardInterrupt:
                mimick_keyboardterminate()
    n_fq = len(fasta) * 2
    # gzip multiprocessing
    allfastq = glob.glob(os.path.abspath(c.OUT) + '/*.fastq')
    chunk_size = len(allfastq)/c.threads
    slices = Chunks(allfastq, math.ceil(chunk_size))
    processes = []

    for i,sli in enumerate(slices):
        #for _s in sli:
        #    mimick_console.log(f'Compressing [blue]{os.path.basename(_s)}[/]')
        p = multiprocessing.Process(target=BGzipper, args=(sli,))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
        
    c.barcodelist.close()
    mimick_console.log(f'Done\n')

if __name__ =='__main__':
    try:
        mimick()
    except KeyboardInterrupt:
        mimick_keyboardterminate()
