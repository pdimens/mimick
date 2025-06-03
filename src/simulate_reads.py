#! /usr/bin/env python3

import os
import sys
import re
import math
import multiprocessing
import time
from tempfile import mkstemp
from .file_ops import *
from .classes import *
from .common import *
from .fastq_writer import *
from .long_molecule import *
import pysam
from pywgsim import wgsim

def wrap_wgsim(*args, **kwargs):
    '''prevents wgsim.core output from being printed to the terminal'''
    # Create temp files for stdout and stderr
    stdout_fd, stdout_path = mkstemp()
    # Save original file descriptors
    orig_stdout = os.dup(1)
    orig_stderr = os.dup(2)
    try:
        # Redirect stdout/stderr at OS level
        os.dup2(stdout_fd, 1)
        os.dup2(stdout_fd, 2)
        # Call the function
        wgsim.core(*args, **kwargs)
    finally:
        # Restore original stdout/stderr
        os.dup2(orig_stdout, 1)
        os.dup2(orig_stderr, 2)
        # Close the temp files
        os.close(stdout_fd)
        os.close(orig_stdout)
        os.close(orig_stderr)
    # Read the captured output
    with open(stdout_path, 'r') as f:
        stdout_content = f.read()
    # Clean up temp files
    os.unlink(stdout_path)
    return stdout_content

def linked_simulation(wgsimparams: wgsimParams,long_molecule: LongMoleculeRecipe, haplotype: int) -> None:
    '''
    The real heavy-lifting that uses pywgsim to simulate short reads from a long molecule that was created and stored
    as a LongMoleculeRecipe. The output FASTQ files are sent to a separate worker thread to append to the final FASTQ
    files without incurring a data race.
    '''
    R1tmp = os.path.abspath(f"{wgsimparams.outdir}/hap{haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R1.tmp")
    R2tmp = os.path.abspath(f"{wgsimparams.outdir}/hap{haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R2.tmp")
    try:
        stdout = wrap_wgsim(
            r1 = R1tmp,
            r2 =R2tmp,
            ref = long_molecule.fasta,
            err_rate = wgsimparams.error,
            mut_rate = wgsimparams.mutation,
            indel_frac = wgsimparams.indels,
            indel_ext = wgsimparams.extindels,
            N = long_molecule.read_count,
            dist = wgsimparams.read_distance,
            stdev = wgsimparams.distance_stdev,
            size_l = wgsimparams.length_R1,
            size_r = wgsimparams.length_R2,
            max_n = 0.05,
            is_hap = 0,
            is_fixed = 0,
            seed = wgsimparams.randomseed
        )
        with open(f'{wgsimparams.outdir}/{wgsimparams.prefix}.wgsim.mutations', 'a') as f:
            f.write(stdout)
    except KeyboardInterrupt:
        mimick_keyboardterminate()
    finally:
        os.remove(long_molecule.fasta)

    if os.stat(R1tmp).st_size == 0:
        os.remove(R1tmp)
        os.remove(R2tmp)
    else:
        WRITER_QUEUE.put((R1tmp, R2tmp, long_molecule.barcode, long_molecule.output_barcode))

