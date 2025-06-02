#! /usr/bin/env python3

import os
import sys
import re
import math
import multiprocessing
import time
from tempfile import mkstemp
from .classes import *
from .file_ops import *
from .long_molecule import *
from .common import *
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

def linked_simulation(schema: Schema, wgsimparams: wgsimParams,long_molecule: LongMoleculeRecipe, molecular_coverage: int|float, outformat: str, processor: int) -> None:
    '''
    The real heavy-lifting that uses pywgsim to simulate short reads from a long molecule that was created and stored
    as a LongMoleculeRecipe. The output FASTQ files are reformatted to include the barcode stored in `long_molecule` in the
    format given by `outformat`.
    '''
    R1tmp = os.path.abspath(f"{container.OUT}/p{processor}.{long_molecule.mol_id}.{long_molecule.barcode}.R1.tmp.fastq")
    R2tmp = os.path.abspath(f"{container.OUT}/p{processor}.{long_molecule.mol_id}.{long_molecule.barcode}.R2.tmp.fastq")
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
        with open(f'{container.OUT}/{container.PREFIX}.p{processor}.wgsim.mutations', 'a') as f:
            f.write(stdout)
    except KeyboardInterrupt:
        mimick_keyboardterminate()

    os.remove(molecule_fasta)
    if os.stat(R1tmp).st_size == 0:
        os.remove(R1tmp)
        os.remove(R2tmp)
    else:
        R1A = os.path.abspath(f'{container.OUT}/{container.PREFIX}.hap_{container.hapnumber.zfill(3)}.R1.fq')
        R2A = os.path.abspath(f'{container.OUT}/{container.PREFIX}.hap_{container.hapnumber.zfill(3)}.R2.fq')

        with open(R1tmp,'r') as infile, open(R1A,'a') as outfile:
            for name,seq,qual in readfq(infile):
                read = format_linkedread(name, barcodestring, seq, qual, True)
                outfile.write('\n'.join(read) + '\n')
        os.remove(R1tmp)
        with open(R2tmp,'r') as infile, open(R2A,'a') as outfile:
            for name,seq,qual in readfq(infile):
                if outformat == "10x":
                    read = [f'@{name}',seq,'+',qual]
                else:
                    read = format_linkedread(name, barcodestring, seq, qual, False)
                outfile.write('\n'.join(read) + '\n')
        os.remove(R2tmp)


