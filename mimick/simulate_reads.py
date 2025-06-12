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
        wgsim.core(
            long_molecule.fasta,
            r1 = R1,
            r2 = R2,
            err_rate = wgsimparams.error,
            mut_rate = wgsimparams.mutation,
            indel_frac = wgsimparams.indels,
            indel_ext = wgsimparams.extindels,
            max_n = 0.05,
            is_hap = 0,
            N = long_molecule.read_count,
            dist = wgsimparams.read_distance,
            stdev = wgsimparams.distance_stdev,
            size_l = wgsimparams.length_R1,
            size_r = wgsimparams.length_R2,
            seed = wgsimparams.randomseed,
            gff = gff
        )

    except KeyboardInterrupt:
        mimick_keyboardterminate()
    except Exception as e:
        return f"{e}"
    finally:
        os.remove(long_molecule.fasta)

    if os.stat(R1).st_size == 0:
        os.remove(R1)
        os.remove(R2)
        os.remove(gff)
        return 0

    append_queue.put((R1, R2, gff, long_molecule.barcode, long_molecule.output_barcode))
    return schema.haplotype, long_molecule.read_count
