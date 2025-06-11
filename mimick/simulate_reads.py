#! /usr/bin/env python3

import os
import re
import sys
from pywgsim import wgsim
from .file_ops import *
from .classes import *
from .common import *
from .process_outputs import *
from .long_molecule import *

def call_with_redirect(func, *args, stdout_file=None, stderr_file=None, **kwargs):
    # Save original file descriptors
    original_stdout = os.dup(1)
    original_stderr = os.dup(2)
    
    try:
        if stdout_file:
            stdout_fd = os.open(stdout_file, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(stdout_fd, 1)
            os.close(stdout_fd)
        
        if stderr_file:
            stderr_fd = os.open(stderr_file, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(stderr_fd, 2)
            os.close(stderr_fd)
        
        # Call the function
        result = func(*args, **kwargs)
        return result
        
    finally:
        # Restore original file descriptors
        os.dup2(original_stdout, 1)
        os.dup2(original_stderr, 2)
        os.close(original_stdout)
        os.close(original_stderr)

def process_recipe(long_molecule: LongMoleculeRecipe, schema: Schema, outdir: str, prefix: str):
    """
    Processess the LongMoleculeRecipe and creates the resulting FASTA file from it using the file connection
    to fasta (the argument). Returns `long_molecule` with the number of reads and the name of the fasta created
    added to it.
    """
    fasta_seq = schema.sequence[max(0,long_molecule.start-1):long_molecule.end+1]
    fasta_header = f'>HAP:{schema.haplotype}_CHROM:{long_molecule.chrom}_START:{long_molecule.start}_END:{long_molecule.end}_BARCODE:{long_molecule.barcode}'
    tempdir = os.path.join(outdir, "temp")
    molecule_fasta = f'{tempdir}/{prefix}_{long_molecule.barcode}.{long_molecule.mol_id}.fa'
    with open(molecule_fasta, 'w') as faout:
        faout.write(
            "\n".join([fasta_header, '\n'.join(re.findall('.{1,60}', fasta_seq))]) + "\n"
        )
    long_molecule.fasta = molecule_fasta
    return long_molecule


def linked_simulation(wgsimparams: wgsimParams, schema: Schema, long_molecule: LongMoleculeRecipe, append_queue) -> None:
    '''
    The real heavy-lifting that uses pywgsim to simulate short reads from a long molecule that was created and stored
    as a LongMoleculeRecipe. The output FASTQ files are sent to a separate worker thread to append to the final FASTQ
    files without incurring a data race.
    '''
    long_molecule = process_recipe(long_molecule, schema, wgsimparams.outdir, wgsimparams.prefix)
    tempdir = os.path.join(wgsimparams.outdir,"temp")
    R1 = f"{tempdir}/hap{schema.haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R1"
    R2 = f"{tempdir}/hap{schema.haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.R2"
    gff = f"{tempdir}/hap{schema.haplotype}.{long_molecule.mol_id}.{long_molecule.barcode}.gff"
    try:
        call_with_redirect(
            wgsim.core, 
            long_molecule.fasta,
            R1,
            R2,
            wgsimparams.error,
            wgsimparams.mutation,
            wgsimparams.indels,
            wgsimparams.extindels,
            0.05,
            0,
            long_molecule.read_count,
            wgsimparams.read_distance,
            wgsimparams.distance_stdev,
            wgsimparams.length_R1,
            wgsimparams.length_R2,
            wgsimparams.randomseed,
            stdout_file=gff,
            stderr_file=os.devnull
        )

    except KeyboardInterrupt:
        mimick_keyboardterminate()
    except Exception as e:
        f"{e}"
    finally:
        os.remove(long_molecule.fasta)

    if os.stat(R1).st_size == 0:
        os.remove(R1)
        os.remove(R2)
        os.remove(gff)
        return 0

    append_queue.put((R1, R2, gff, long_molecule.barcode, long_molecule.output_barcode))
    return schema.haplotype, long_molecule.read_count
