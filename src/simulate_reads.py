#! /usr/bin/env python3

import os
import re
import sys
import subprocess
from .file_ops import *
from .classes import *
from .common import *
from .process_outputs import *
from .long_molecule import *
import pysam

def process_recipe(fasta: pysam.FastaFile, long_molecule: LongMoleculeRecipe, schema: Schema, rng):
    """
    Processess the LongMoleculeRecipe and creates the resulting FASTA file from it using the file connection
    to fasta (the argument). Returns `long_molecule` with the number of reads and the name of the fasta created
    added to it.
    """
    fasta_seq = schema.sequence[long_molecule.start-1:long_molecule.end+1]
    fasta_header = f'>CHROM:{long_molecule.chrom}_START:{long_molecule.start}_END:{long_molecule.end}_BARCODE:{long_molecule.barcode}'
    molecule_fasta = os.path.abspath(f'{long_molecule.out_prefix}_{long_molecule.barcode}.{long_molecule.mol_id}.fa')
    with open(molecule_fasta, 'w') as faout:
        faout.write(
            "\n".join([fasta_header, '\n'.join(re.findall('.{1,60}', fasta_seq))]) + "\n"
        )

    normalized_length = len(fasta_seq)-fasta_seq.count('N')
    if schema.mol_coverage < 1:
        N = max(1,int(normalized_length * schema.mol_coverage)/(schema.read_length*2))
    else:
        # draw N from a normal distribution with a mean of molcov and stdev of molcov/3, avoiding < 0
        N = max(0, rng.normal(schema.mol_coverage, schema.mol_coverage/3))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, normalized_length/(schema.read_length*2))
    if schema.singletons > 0 and N != 0:
        if rng.uniform(0,1) > schema.singletons:
            N = 1
    long_molecule.fasta = molecule_fasta
    long_molecule.read_count = int(N)
    return long_molecule

def linked_simulation(fasta, wgsimparams: wgsimParams, schema: Schema, long_molecule: LongMoleculeRecipe, rng, append_queue) -> None:
    '''
    The real heavy-lifting that uses pywgsim to simulate short reads from a long molecule that was created and stored
    as a LongMoleculeRecipe. The output FASTQ files are sent to a separate worker thread to append to the final FASTQ
    files without incurring a data race.
    '''
    long_molecule = process_recipe(fasta, long_molecule, schema, rng)
    R1 = os.path.abspath(f"{wgsimparams.outdir}/hap{schema.haplotype_number}.{long_molecule.mol_id}.{long_molecule.barcode}.R1")
    R2 = os.path.abspath(f"{wgsimparams.outdir}/hap{schema.haplotype_number}.{long_molecule.mol_id}.{long_molecule.barcode}.R2")
    gff = os.path.abspath(f"{wgsimparams.outdir}/hap{schema.haplotype_number}.{long_molecule.mol_id}.{long_molecule.barcode}.gff")
    try:
        with open(gff, "w") as _gff:
            subprocess.check_call(
                [
                    "pywgsim",
                    "--err",
                    str(wgsimparams.error),
                    "--mut",
                    str(wgsimparams.mutation),
                    "--frac",
                    str(wgsimparams.indels),
                    "--ext",
                    str(wgsimparams.extindels),
                    "-N",
                    str(long_molecule.read_count),
                    "--dist",
                    str(wgsimparams.read_distance),
                    "--stdev",
                    str(wgsimparams.distance_stdev),
                    "--L1",
                    str(wgsimparams.length_R1),
                    "--L2",
                    str(wgsimparams.length_R2),
                    "--amb",
                    "0.05",
                    "--seed",
                    str(wgsimparams.randomseed),
                    long_molecule.fasta,
                    R1,
                    R2
                ],
                stdout= _gff
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
    else:
        append_queue.put((R1, R2, gff, long_molecule.barcode, long_molecule.output_barcode))
        return long_molecule.read_count
