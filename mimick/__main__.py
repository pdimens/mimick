#!/usr/bin/env python3

import os
import sys
import gzip
import threading
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from random import choices, getrandbits
from time import sleep
import numpy as np
import rich_click as click
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
            "options": ["--help", "--circular", "--output-prefix", "--output-type", "--quiet", "--regions", "--seed", "--threads", "--version"],
            "panel_styles": {"border_style": "dim"}
        },
        {
            "name": "Read Simulation Parameters",
            "options": ["--coverage","--distance","--error","--extindels","--indels","--lengths","--mutation","--stdev"],
            "panel_styles": {"border_style": "dim blue"}
        },
        {
            "name": "Linked Read Parameters",
            "options": ["--molecule-attempts", "--molecule-coverage", "--molecule-length", "--molecules-per", "--segments", "--singletons"],
            "panel_styles": {"border_style": "dim magenta"}
        },
    ]
}

@click.version_option("0.0.0", prog_name="mimick")
@click.command(epilog = "Documentation: https://pdimens.github.io/mimick/", no_args_is_help = True)
@click.option('-C','--circular', is_flag= True, default = False, help = 'contigs are circular/prokaryotic')
@click.option('-o','--output-prefix', help='output file prefix', type = click.Path(exists = False, writable=True, resolve_path=True), default = "simulated/SIM", show_default=True)
@click.option('-O','--output-type', help='FASTQ output format', type = click.Choice(["10x", "stlfr", "standard", "standard:haplotagging", "standard:stlfr", "haplotagging", "tellseq"], case_sensitive=False))
@click.option('-q','--quiet', show_default = True, default = "0", type = click.Choice([0,1,2]), help = '`0` all output, `1` no progress bar, `2` no output')
@click.option('-r','--regions', help='one or more regions to simulate, in BED format', type = click.Path(dir_okay=False, readable=True, resolve_path=True))
@click.option('-t','--threads', help='number of threads to use for simulation', type=click.IntRange(min=1), default=2, show_default=True)
@click.option('-S','--seed', help='random seed for simulation', type=click.IntRange(min=1), default=None)
#Paired-end FASTQ simulation using pywgsim
@click.option('--coverage', help='mean coverage (depth) target for simulated data', show_default=True, default=30.0, type=click.FloatRange(min=0.05))
@click.option('--distance', help='outer distance between the two read ends in bp', default=500, show_default=True, type=click.IntRange(min=0))
@click.option('--error', help='base error rate', default=0.02, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--extindels', help='indels extension rate', default=0.25, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--indels', help='indels rate', default=0.15, show_default=True, type=click.FloatRange(min=0, max=1))
@click.option('--lengths', help='length of R1,R2 reads in bp', default="150,150", show_default=True, type=ReadLengths())
@click.option('--mutation', help='mutation rate', default=0.001, show_default=True, type=click.FloatRange(min=0))
@click.option('--stdev', help='standard deviation for `--distance`', default=50, show_default=True, type=click.IntRange(min=0))
#Linked-read simulation
@click.option('-a','--molecule-attempts', help='how many tries to create a molecule with <70% ambiguous bases', show_default=True, default=300, type=click.IntRange(min=5))
@click.option('-c','--molecule-coverage', help='mean percent coverage per molecule if <1, else mean number of reads per molecule', default=0.2, show_default=True, type=click.FloatRange(min=0.00001))
@click.option('-m','--molecule-length', help='mean length of molecules in bp', show_default=True, default=80000, type=click.IntRange(min=650))
@click.option('-n','--molecules-per', help='mean number of unrelated molecules per barcode per chromosome, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Normal distribution', default=2, show_default=True, type=int)
@click.option('-x', '--segments', help='treat barcodes as combinatorial with this many segments', default = 4, show_default=True, type= click.IntRange(min=1))
@click.option('-s','--singletons', help='proportion of barcodes that will only have one read pair', default=0, show_default=True, type=click.FloatRange(0,1))
@click.argument('barcodes', type = Barcodes())
@click.argument('fasta', type = click.Path(exists=True, dir_okay=False, resolve_path=True, readable=True), nargs = -1, required=True)
def mimick(barcodes, fasta, circular, output_prefix, output_type, quiet, seed, regions, threads,coverage,distance,error,extindels,indels,lengths,mutation,stdev,segments, molecule_coverage, molecule_attempts, molecule_length, molecules_per, singletons):
    """
    Simulate linked-read FASTQ using genome haplotypes. Barcodes can be supplied one of two ways:

    1. let Mimick randomly generate barcodes based on a specification of `length,count`
        - two integers, comma-separated, no space
        - e.g. `16,400000` would generate 400,000 unique 16bp barcodes 
    2. you can provide a file of specific nucleotide barcodes, 1 per line

    In addition to selecting an `--output-type`, barcodes can be parsed absolutely or you can specify the
    linked-read barcode type using `-x/--segments`. For example, to simulate the common 4-segment haplotagging style,
    use `-x 4` and have `--output-type haplotagging` (or use a different output style, if preferred). The `standard` output
    types can be suffixed with `:haplotagging` or `:stlfr` to use those barcode styles with the standard format
    (e.g. `standard:haplotagging`). The table below serves as a guide for the configurations for the common linked-read varieties: 

    | chemistry    | `--segments` | `--lengths` | Format                                       |
    |:-------------|:------------:|:-----------:|:---------------------------------------------|
    | 10x/tellseq  |     `1`      |  `134,150`  | single barcode on R1                         |
    | haplotagging |     `4`      |  `150,150`  | I1 and I2 each with combinatorial 2-barcodes |
    | stlfr        |     `3`      |  `150,108`  | combinatorial 3-barcode on R2                |

    | --output-type    | Barcode Location              | default for | Example                    |
    |:-----------------|:------------------------------|:-----------:|:---------------------------|
    | `10x`            | start of R1 sequence          |             | `ATAGACCATAGA`GGACA...     |
    | `haplotagging`   | header comment as `BX:Z:ACBD` |             | `@SEQID BX:Z:A0C331B34D87` |
    | `standard[:...]` | `BX:Z:BARCODE VX:i:N`         | all others  | `@SEQID BX:Z:ATACGAGACA`   |
    | `stlfr`          | sequence ID + `#1_2_3`        |   `-x 3`    | `@SEQID#1_354_39`          |
    | `tellseq`        | sequence ID + `:ATCG`         |   `-x 1`    | `@SEQID:TATTAGCAC`         |
    """
    if molecules_per == 0:
        error_terminate("The value for [yellow]--molecule-number[/] cannot be 0.")
    if regions and circular:
        error_terminate(
            "[yellow]--circular[/] cannot be used with [yellow]--regions[/] beacuse Mimick may incorrectly circularize"
            " contig slices. If you are using [yellow]--regions[/] to simulate from a subset of contigs, but in their"
            " entirety, then create new FASTA files containing only those contigs of interest and use [yellow]--circular[/]"
            " without [yellow]--regions[/]."
        )
    PROGRESS.disable = quiet > 0
    BARCODE_OUTPUT_FORMAT = output_type.lower() if output_type else OUTTYPE.get(segments, "standard")

    RNG = np.random.default_rng(seed = seed)
    WGSIMPARAMS = wgsimParams(
        error,
        mutation,
        indels,
        extindels,
        distance,
        stdev,
        lengths[0],
        lengths[1],
        seed if seed else getrandbits(16),
        output_prefix
    )
    tempdir = os.path.join(WGSIMPARAMS.outdir, "temp")
    os.makedirs(os.path.join(tempdir, "molecules"), exist_ok = True)

    if isinstance(barcodes, str):
        BARCODE_PATH = barcodes
        if quiet < 2:
            mimick_console.log(f'Validating barcodes in [blue]{os.path.basename(BARCODE_PATH)}[/]')
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
    BARCODES = interpret_barcodes(BARCODE_PATH, segments, BARCODE_OUTPUT_FORMAT)
    BARCODES.remaining = BARCODES.max

    if quiet == 0:
        PROGRESS.start()
    # create a countdown progressbar tracking remaining barcodes
    _progress_bc = PROGRESS.add_task(f"[bold]Barcodes Remaining", total=BARCODES.remaining, completed= BARCODES.remaining)

    _index_fasta = PROGRESS.add_task(f"[bold magenta]Process inputs", total=len(fasta))

    fasta_indexes = []
    for idx,haplotype in enumerate(fasta):
        if quiet < 2:
            mimick_console.log(f"Indexing [{STYLES[idx]}]{os.path.basename(haplotype)}[/]")
        fasta_indexes += index_fasta(haplotype)
        PROGRESS.update(_index_fasta, advance=1)
    PROGRESS.update(_index_fasta, visible=False)
    # to be more accurate with target read calculation, get the average read length after accounting for the sequencing
    # space occupied by the barcodes
    avg_adjusted_read_length = (WGSIMPARAMS.length_R1 + WGSIMPARAMS.length_R2) / 2
    if regions:
        if quiet < 2:
            mimick_console.log(f'Processing haplotypes and intervals')
        SIMULATION_SCHEMA = BEDtoInventory(regions, fasta, coverage/len(fasta), molecule_coverage, molecule_length, avg_adjusted_read_length, molecule_length, singletons, circular)
    else:
        if quiet < 2:
            mimick_console.log(f'Processing haplotypes')
        SIMULATION_SCHEMA = FASTAtoInventory(fasta, coverage/len(fasta), molecule_coverage, molecule_length, avg_adjusted_read_length, singletons, circular)

    # indices are no longer needed, so remove them
    for fai in fasta_indexes:
        os.remove(fai)

    # initialize total progress bar
    _progress_sim = PROGRESS.add_task(f"[blue]Total Progress", total=sum([i.reads_required for i in SIMULATION_SCHEMA.values()]))

    # initialize progress bar for each haplotype
    _progress_haplo = []
    for i,fa in enumerate(fasta,1):
        _progress_haplo.append(PROGRESS.add_task(f"[{STYLES[i-1]}]Haplotype {i}", total=sum([j.reads_required for j in SIMULATION_SCHEMA.values() if j.haplotype == i])))

    # this list will be iterated over and have entries removed when their target read counts are met
    schemas = list(SIMULATION_SCHEMA.keys())
    quota_reached = set()

    # housekeeping for multithreading
    max_queued = max(10000, 1000 * threads)
    futures = set()

    # Create a thread-safe queue to write temp files to final outputs
    # setup the initialization arguments (final output file names, output type, verbosity, etc.)
    output_appender = FileProcessor(output_prefix, BARCODE_OUTPUT_FORMAT, quiet)
    output_appender.start()
    if quiet < 2:
        mimick_console.log(f'Simulating molecules and reads')

    # file to record all the molecules that were created
    MOLECULE_INVENTORY = MoleculeRecorder(output_prefix)

    with ThreadPoolExecutor(max_workers = threads - 1) as executor:
        # number of molecules is fixed if option was < 0
        n_molecules = abs(molecules_per)
        # only need to calculate the stdev once in the event n_molecules needs to be sampled in the main loop
        sd = molecules_per/(molecules_per - 2) if molecules_per > 2 else 3/4
        # finally, iterate over the barcode generator until all quotas in the schema are met
        # the while loop breaks when there are no more schema left (schema get removed when their targets are achieved)
        while True:
            try:
                selected_bc = BARCODES.get_next_bc()
                BARCODES.remaining -= 1
                PROGRESS.update(_progress_bc, completed= BARCODES.remaining)
            except StopIteration:
                stop_event.set()
                executor.shutdown(wait = False, cancel_futures=True)
                output_appender.stop()
                PROGRESS.stop()
                error_terminate('No more barcodes left for simulation. The requested parameters require more barcodes.')

            output_bc = BARCODES.get_next_out(selected_bc)

            # if necessary, draw n molecules from a normal distribution
            if molecules_per > 0:
                n_molecules = max(1, int(RNG.normal(molecules_per, sd)))

            # randomly pick a number of intervals equal to the number of molecules associated with the barcode
            # sampled with replacement, so e.g. it's possible that when n_molecules = 3, 2 of the unrelated molecules come from
            # interval 1, and 1 molecule comes from interval 2. In other words, this method (crucially!) randomizes
            # what interval/chromosome the unrelated molecules are drawn from instead of the molecules always being from the same
            # chromosome/interval
            for target in choices(schemas, k = n_molecules):
                _haplotype = SIMULATION_SCHEMA[target].haplotype
                molecule_recipe = create_long_molecule(SIMULATION_SCHEMA[target], RNG, selected_bc, output_bc, WGSIMPARAMS, molecule_attempts)
                if not molecule_recipe:
                    executor.shutdown(wait = False, cancel_futures=True)
                    error_terminate((
                        f"Unable to create a molecule with <70% ambiguous ([bold yellow]N[/]) bases after {molecule_attempts} tries for [default]"
                        f"{SIMULATION_SCHEMA[target].chrom}:{SIMULATION_SCHEMA[target].start}:{SIMULATION_SCHEMA[target].end} (haplotype {SIMULATION_SCHEMA[target].haplotype})[/]."
                        f" Consider increasing [blue]--molecule-attempts[/] or providing intervals to [blue]--regions[/] that avoid very long stretches of ambiguous "
                        "bases, which is likely the cause of this error."), appender=output_appender)
                
                MOLECULE_INVENTORY.write(molecule_recipe)
                # add number of reads to the tracker in the schema
                    
                future = executor.submit(linked_simulation, WGSIMPARAMS, molecule_recipe)
                futures.add(future)

                SIMULATION_SCHEMA[target].reads_current += molecule_recipe.read_count
                # if the number of simulated reads reaches the target, add the index to the "remove from schema" list
                if SIMULATION_SCHEMA[target].reads_current >= SIMULATION_SCHEMA[target].reads_required:
                    quota_reached.add(target)
            
            # check to see if read targets were met and if so, remove them from consideration in the next iteration
            if quota_reached:
                schemas = list(set(schemas) - quota_reached)
                quota_reached = set()

            # update progress with backpressure control
            while True:
                done = [f for f in futures if f.done()]
                for f in done:
                    futures.remove(f)
                    _result = f.result()
                    if _result== 0:
                        continue
                    if isinstance(_result, str):
                        executor.shutdown(wait = False, cancel_futures=True)
                        error_terminate(_result, appender=output_appender)
                    _hap = _result.haplotype
                    _N = _result.read_count
                    output_appender.submit_files(_result)
                    PROGRESS.update(_progress_sim, advance =_N)
                    PROGRESS.update(_progress_haplo[_hap-1], advance = _N)
                if len(futures) < max_queued:
                    break
                sleep(0.05)

            # if the target coverage was reached for all intervals, we can move on to the next haplotype
            if not schemas:
                break

        # finish processing any molecules that were in the thread queue
        for f in as_completed(futures):
            _result = f.result()
            if _result == 0:
                continue
            if isinstance(_result, str):
                executor.shutdown(wait = False, cancel_futures=True)
                error_terminate(_result, appender=output_appender)
            output_appender.submit_files(_result)
            _hap = _result.haplotype
            _N = _result.read_count
            PROGRESS.update(_progress_sim, advance = _N)
            PROGRESS.update(_progress_haplo[_hap-1], advance = _N)
        # wait for the final-output-writer thread to finish
        output_appender.stop()
    shutil.rmtree(tempdir, ignore_errors=True)
    MOLECULE_INVENTORY.close()
    PROGRESS.stop()
    if quiet < 2:
        mimick_console.log("Finished!")
        #mimick_console.rule("[dim]Separate by haplotype[/]", style="dim")
        #mimick_console.print("[dim]If you would like to separate the reads by haplotype, use:[/]")
        #mimick_console.print(f"    [dim green]zgrep -A3 \"^@HAP:X_\" {os.path.relpath(output_prefix)}.R1.fq.gz[/]")
        #mimick_console.print("[dim]Where [dim green]X[/] is the haplotype you are trying to isolate.[/]")

if __name__ =='__main__':
    try:
        mimick()
    except KeyboardInterrupt:
        mimick_keyboardterminate()
    except Exception as e:
        error_terminate(e)
        PROGRESS.stop()
