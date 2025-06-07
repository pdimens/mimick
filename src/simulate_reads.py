#! /usr/bin/env python3

import os
import sys
from contextlib import contextmanager
from .file_ops import *
from .classes import *
from .common import *
from .process_outputs import *
from .long_molecule import *
import pysam
from pywgsim import wgsim

@contextmanager
def redirect_stdout_stderr(stdout_path, stderr_path):
    # Save original file descriptors
    original_stdout_fd = sys.stdout.fileno()
    original_stderr_fd = sys.stderr.fileno()

    # Open the output files
    with open(stdout_path, 'w') as out, open(stderr_path, 'w') as err:
        # Flush any existing content in stdout/stderr
        sys.stdout.flush()
        sys.stderr.flush()

        # Save copies of original fds
        saved_stdout_fd = os.dup(original_stdout_fd)
        saved_stderr_fd = os.dup(original_stderr_fd)

        try:
            # Redirect stdout and stderr to the file
            os.dup2(out.fileno(), original_stdout_fd)
            os.dup2(err.fileno(), original_stderr_fd)
            yield
        finally:
            # Restore original fds
            os.dup2(saved_stdout_fd, original_stdout_fd)
            os.dup2(saved_stderr_fd, original_stderr_fd)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)

def linked_simulation(wgsimparams: wgsimParams,long_molecule: LongMoleculeRecipe, haplotype: int, append_queue) -> None:
    '''
    The real heavy-lifting that uses pywgsim to simulate short reads from a long molecule that was created and stored
    as a LongMoleculeRecipe. The output FASTQ files are sent to a separate worker thread to append to the final FASTQ
    files without incurring a data race.
    '''
    R1 = os.path.abspath(f"{wgsimparams.outdir}/hap{haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R1.tmp")
    R2 = os.path.abspath(f"{wgsimparams.outdir}/hap{haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R2.tmp")
    gff = os.path.abspath(f"{wgsimparams.outdir}/hap{haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.gff.tmp")
    try:
        stdout_file = gff
        stderr_file = f"{gff}.stderr"
        with redirect_stdout_stderr(stdout_file, stderr_file):
            wgsim.core(
                r1 = R1,
                r2 =R2,
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
    except KeyboardInterrupt:
        mimick_keyboardterminate()
    finally:
        os.remove(long_molecule.fasta)
        os.remove(stderr_file)

    if os.stat(R1).st_size == 0:
        os.remove(R1)
        os.remove(R2)
    else:
        append_queue.put((R1, R2, gff, long_molecule.barcode, long_molecule.output_barcode))
        return long_molecule.read_count
