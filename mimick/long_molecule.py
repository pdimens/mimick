#! /usr/bin/env python3

import os
from random import getrandbits
import numpy as np
from .classes import Schema

class LongMoleculeRecipe(object):
    '''Molecule instance'''
    def __init__(self, haplotype, fasta, chrom, start, end, barcode,outbarcode,mol_id, read_count):
        self.haplotype = haplotype
        self.fasta = fasta
        self.output_basename = f"hap{haplotype}.{mol_id}.{barcode}"
        self.chrom=chrom
        self.start=start
        self.end=end
        self.barcode=barcode
        self.output_barcode=outbarcode
        self.mol_id=mol_id
        self.read_count = int(read_count)
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            outstring += f"{i}: {j}\n"
        return outstring

def create_long_molecule(schema: Schema, rng, barcode, outputbarcode, wgsimparams) -> LongMoleculeRecipe:
    '''
    Randomly generates a long molecule and writes it to a FASTA file.
    Length of molecules is randomly distributed using an exponential distribution, with a minimum of 650bp.
    Returns a LongMoleculeRecipe that contains all the necessary information to simulate reads from that molecule.
    '''
    molnumber = getrandbits(32)
    len_interval = schema.end - schema.start
    # make sure to cap the molecule length to the length of the interval/chromosome
    molecule_length = 0
    # make sure the molecule is greater than 650
    while molecule_length < 650 or molecule_length > len_interval:
        molecule_length = rng.exponential(scale = schema.mol_length)

    # set the max position to be length - mol_length to avoid additional computation
    molecule_length = int(molecule_length)
    start = int(rng.uniform(low = 0, high = len_interval - molecule_length))
    end = start + molecule_length - 1

    fasta_seq = schema.sequence[start:end+1]
    normalized_length = len(fasta_seq)-fasta_seq.count('N')

    # if a singleton proportion is provided, conditionally drop the number of reads to 1
    if schema.singletons > 0 and rng.uniform(0,1) > schema.singletons:
        N = 1
    elif schema.mol_coverage < 1:
        # set a minimum number of reads to 2 to avoid singletons
        N = max(2, normalized_length * schema.mol_coverage)/(schema.read_length*2)
    else:
        # draw N from an exponential distribution with a minimum set to 2 reads to avoid singletons
        N = max(2, rng.exponential(schema.mol_coverage))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, normalized_length/(schema.read_length*2))

    tempdir = os.path.join(wgsimparams.outdir, "temp","molecules")
    fasta_file = f'{tempdir}/{wgsimparams.prefix}_{barcode}.{molnumber}.fa'
    fasta_header = f'>HAP:{schema.haplotype}_CHROM:{schema.chrom}_START:{start}_END:{end}_BARCODE:{barcode}'

    with open(fasta_file, 'w') as faout:
        faout.write(
            "\n".join([fasta_header, fasta_seq])
        )
    return LongMoleculeRecipe(schema.haplotype, fasta_file, schema.chrom, start, end, barcode, outputbarcode, molnumber, N)
