#!/usr/bin/env python3

import glob
import os
import sys
import re
import math
import gzip
import multiprocessing
import threading
from itertools import product
from concurrent.futures import ThreadPoolExecutor, as_completed
from random import choices, getrandbits
from time import sleep
import numpy as np
import rich_click as click
from rich.progress import Progress, TextColumn, TimeElapsedColumn, TaskProgressColumn, BarColumn
import pysam
from .classes import *
from .cli_classes import *
from .common import *
from .file_ops import *
from .long_molecule import *
from .simulate_reads import *

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

PROGRESS = Progress(
    TextColumn("[progress.description]{task.description}"),
    BarColumn(finished_style="purple", complete_style="yellow"),
    TaskProgressColumn(),
    TimeElapsedColumn(),
    transient=True,
    console=mimick_console
)
stop_event = threading.Event()

@click.version_option("0.0.0", prog_name="mimick")
@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/SIM", show_default=True)
@click.option('-O','--output-type', help='output format of FASTQ files', type = click.Choice(["10x", "stlfr", "standard", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-q','--quiet', show_default = True, default = "0", type = click.Choice([0,1,2]), help = '`0` all output, `1` no progress bar, `2` no output')
@click.option('-r','--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('-t','--threads', help='number of threads to use for simulation', type=click.IntRange(min=1), default=2, show_default=True)
@click.option('-S','--seed', help='random seed for simulation', type=click.IntRange(min=1), default=None)
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
    if molecules_per == 0:
        error_terminate("The value for [yellow]--molecule-number[/] cannot be 0.")

    PROGRESS.disable = quiet > 0
    LR_CHEMISTRY = lr_type.lower()
    BARCODE_OUTPUT_FORMAT = output_type.lower() if output_type else LR_CHEMISTRY
    RNG = np.random.default_rng(seed = seed)
    WGSIMPARAMS = wgsimParams(
        error,
        mutation,
        indels,
        extindels,
        distance,
        stdev,
        0,
        0,
        seed if seed else getrandbits(16),
        output_prefix
    )

    os.makedirs(WGSIMPARAMS.outdir, exist_ok= True)

    if isinstance(barcodes, str):
        BARCODE_PATH = barcodes
        if quiet < 2:
            mimick_console.log(f'Validating barcodes in [blue]{os.path.basename(barcodes)}[/]')
        try:
            with gzip.open(BARCODE_PATH, 'rt') as filein:
                BARCODES, BARCODE_LENGTH_BP, BARCODES_TOTAL_COUNT = interpret_barcodes(filein, LR_CHEMISTRY)
        except gzip.BadGzipFile:
            with open(BARCODE_PATH, 'r') as filein:
                BARCODES, BARCODE_LENGTH_BP, BARCODES_TOTAL_COUNT = interpret_barcodes(filein, LR_CHEMISTRY)
        except:
            error_terminate(f'Cannot open {BARCODE_PATH} for reading')
    else:
        BARCODE_PATH = f"{output_prefix}.generated.barcodes"
        bp,count = barcodes
        if quiet < 2:
            mimick_console.log(f'Generating {count} barcodes at [default bold]{bp}bp[/] each')
        with open(BARCODE_PATH, "w") as bc_out:
            for i,bc in enumerate(generate_barcodes(bp),1):
                bc_out.write("".join(bc) + "\n")
                if i == count:
                    break
        with open(BARCODE_PATH, 'r') as filein:
            BARCODES, BARCODE_LENGTH_BP, BARCODES_TOTAL_COUNT = interpret_barcodes(filein, LR_CHEMISTRY)

    BARCODES_REMAINING = BARCODES_TOTAL_COUNT

    # process the read lengths and modify them to account for the type of chemistry desired
    if LR_CHEMISTRY in ["10x", "tellseq"]:
        # barcode at beginning of read 1
        WGSIMPARAMS.length_R1 = length - BARCODE_LENGTH_BP
        if WGSIMPARAMS.length_R1 <= 5:
            error_terminate(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {length}, Barcode length: {BARCODE_LENGTH_BP}')
        WGSIMPARAMS.length_R2 = length
    elif LR_CHEMISTRY == "stlfr":
        # barcode at the end of read 2
        WGSIMPARAMS.length_R1 = length
        WGSIMPARAMS.length_R2 = length - BARCODE_LENGTH_BP
        if WGSIMPARAMS.length_R2 <= 5:
            error_terminate(f'Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {length}, Barcode length: {BARCODE_LENGTH_BP}')
    else:
        # would be 4-segment haplotagging where AC on read 1 and BD on read 2, but in index read, so read not shortened
        WGSIMPARAMS.length_R1 = length
        WGSIMPARAMS.length_R2 = length

    # setup barcode generator #
    if BARCODE_OUTPUT_FORMAT == "haplotagging":
        if BARCODES_TOTAL_COUNT > 96**4:
            error_terminate(f'The barcodes and barcode type supplied will generate a potential {BARCODES_TOTAL_COUNT} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')
        bc_range = [f"{i}".zfill(2) for i in range(1,97)]
        BARCODE_OUTPUT_GENERATOR = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)

    if BARCODE_OUTPUT_FORMAT == "stlfr":
        if BARCODES_TOTAL_COUNT > 1537**3:
            error_terminate(f'The barcodes and barcode type supplied will generate a potential {BARCODES_TOTAL_COUNT} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')
        bc_range = range(1, 1537)
        BARCODE_OUTPUT_GENERATOR = product(bc_range, bc_range, bc_range)

    # number of molecules is fixed if option was < 0
    n_molecules = abs(molecules_per)
    # only need to calculate the stdev once in the event n_molecules needs to be sampled in the main loop
    sd = molecules_per/(molecules_per - 2) if molecules_per > 2 else 3/4
    
    if quiet == 0:
        PROGRESS.start()
    # create a countdown progressbar tracking remaining barcodes
    _progress_bc = PROGRESS.add_task(f"[bold magenta]Barcodes Remaining", total=BARCODES_REMAINING, completed= BARCODES_REMAINING)
    
    # file to record all the molecules that were created
    MOLECULE_INVENTORY = open(f'{output_prefix}.molecules', 'w')

    # add variation to the console rules. necessary? no. fun? yes.
    styles = ["purple", "yellow", "green", "orange", "blue", "magenta"] * 4

    with ThreadPoolExecutor(max_workers = threads - 1) as executor:
        futures = set()
        max_queued = 2 * threads
        for idx,haplotype in enumerate(fasta, 1):
            if quiet < 2:
                mimick_console.rule(os.path.basename(haplotype), style = styles[idx-1])
            index_fasta(haplotype)
            MOLECULE_INVENTORY.write(
                "\t".join(["haplotype", "chromosome", "start_position", "end_position", "length", "reads", "nucleotide_barcode", "output_barcode"]) + "\n"
            )
            
            if regions:
                if quiet < 2:
                    mimick_console.log(f'Processing sequences and intervals')
                SIMULATION_SCHEMA = BEDtoInventory(regions, haplotype, coverage/len(fasta), molecule_coverage, molecule_length, length, molecule_length, singletons, idx)
            else:
                if quiet < 2:
                    mimick_console.log(f'Processing sequences')
                SIMULATION_SCHEMA = FASTAtoInventory(haplotype, coverage/len(fasta), molecule_coverage, molecule_length, length, singletons, idx)
            
            _progress_sim = PROGRESS.add_task(f"[default]Haplotype {idx}", total=sum([i.reads_required for i in SIMULATION_SCHEMA.values()]))

            # this list will be iterated over and have entries removed when their target read counts are met
            schemas = list(SIMULATION_SCHEMA.keys())
            quota_reached = set()

            # Create a thread-safe queue to write temp files to final outputs
            # setup the initialization arguments (final output file names, output type, verbosity, etc.)
            WRITER_QUEUE = queue.Queue()
            appender_args = (
                os.path.abspath(f'{output_prefix}.hap_{str(idx).zfill(3)}.R1.fq'),
                os.path.abspath(f'{output_prefix}.hap_{str(idx).zfill(3)}.R2.fq'),
                os.path.abspath(f'{output_prefix}.hap_{str(idx).zfill(3)}.gff'),
                BARCODE_OUTPUT_FORMAT,
                quiet,
                WRITER_QUEUE
            )
            output_appender = threading.Thread(target=append_worker, args = appender_args)
            output_appender.start()

            if quiet < 2:
                mimick_console.log(f'Simulating molecules and reads')

            # finally, iterate over the barcode generator until all quotas in the schema are met
            # the while loop breaks when there are no more schema left (schema get removed when their targets are achieved)
            while True:
                try:
                    selected_bc = "".join(next(BARCODES))
                    BARCODES_REMAINING -= 1
                    PROGRESS.update(_progress_bc, completed= BARCODES_REMAINING)
                except StopIteration:
                    stop_event.set()
                    executor.shutdown(wait = False, cancel_futures=True)
                    WRITER_QUEUE.put(None)
                    output_appender.join()
                    PROGRESS.stop()
                    error_terminate('No more barcodes left for simulation. The requested parameters require more barcodes.')

                # setup output barcode formats
                if BARCODE_OUTPUT_FORMAT == "haplotagging":
                    output_bc = "".join(next(BARCODE_OUTPUT_GENERATOR))
                elif BARCODE_OUTPUT_FORMAT == "stlfr":
                    output_bc = "_".join(str(i) for i in next(BARCODE_OUTPUT_GENERATOR))
                else:
                    output_bc = selected_bc

                # if necessary, draw n molecules from a normal distribution
                if molecules_per > 0:
                    n_molecules = max(1, int(RNG.normal(molecules_per, sd)))
                
                # randomly pick a number of intervals equal to the number of molecules associated with the barcode
                # sampled with replacement, so e.g. it's possible that when n_molecules = 3, 2 of the unrelated molecules come from
                # interval 1, and 1 molecule comes from interval 2. In other words, this method (crucially!) randomizes
                # what interval/chromosome the unrelated molecules are drawn from instead of the molecules always being from the same
                # chromosome/interval 
                for target in choices(schemas, k = n_molecules):
                    molecule_recipe = create_long_molecule(SIMULATION_SCHEMA[target], RNG, selected_bc, output_prefix, output_bc)
                    MOLECULE_INVENTORY.write(
                       "\t".join([f"haplotype_{idx}", molecule_recipe.chrom, str(molecule_recipe.start), str(molecule_recipe.end), str(molecule_recipe.end-molecule_recipe.start), str(molecule_recipe.read_count),  selected_bc, output_bc]) + "\n"
                    )
                    SIMULATION_SCHEMA[target].reads_current += molecule_recipe.read_count
                    if SIMULATION_SCHEMA[target].reads_current >= SIMULATION_SCHEMA[target].reads_required:
                        quota_reached.add(target)

                    # Backpressure: wait if too many futures are outstanding
                    while futures and len(futures) >= max_queued:
                        # Poll for completed futures to reduce pressure
                        done = {f for f in futures if f.done()}
                        for f in done:
                            futures.remove(f)
                            PROGRESS.update(_progress_sim, advance=f.result())
                        sleep(.1)
                    future = executor.submit(linked_simulation, WGSIMPARAMS, SIMULATION_SCHEMA[target], molecule_recipe, WRITER_QUEUE)
                    futures.add(future)
                # peek the threads to see if any finished and update the progress bar with any that did
                # this is duplicated to account for checking when we aren't backlogged with submissions
                done = {f for f in futures if f.done()}
                for f in done:
                    futures.remove(f)
                    PROGRESS.update(_progress_sim, advance=f.result())
                # check to see if read targets were met and if so, remove them from consideration in the next iteration
                if quota_reached:
                    schemas = list(set(schemas) - quota_reached)
                    quota_reached = set()
                    # if the target coverage was reached for all intervals, we can move on to the next haplotype        
                    if not schemas:
                        break
            # finish processing any molecules that were in the thread queue
            for f in as_completed(futures):
                PROGRESS.update(_progress_sim, advance=f.result())
            # wait for the final-output-writer thread to finish
            WRITER_QUEUE.put(None)
            output_appender.join()
    MOLECULE_INVENTORY.close()
    PROGRESS.stop()

if __name__ =='__main__':
    try:
        mimick()
    except KeyboardInterrupt:
        mimick_keyboardterminate()
    except Exception as e:
        PROGRESS.stop()
