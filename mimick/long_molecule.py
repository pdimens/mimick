#! /usr/bin/env python3

import os
from random import getrandbits
import numpy as np
from .classes import Schema

class LongMoleculeRecipe(object):
    '''Molecule instance'''
    def __init__(self, chrom, start, end, barcode,outbarcode,mol_id, read_count, out_prefix):
        self.fasta = ""
        self.chrom=chrom
        self.start=start
        self.end=end
        self.barcode=barcode
        self.output_barcode=outbarcode
        self.mol_id=mol_id
        self.out_prefix = out_prefix
        self.read_count = read_count
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            outstring += f"{i}: {j}\n"
        return outstring

def create_long_molecule(schema: Schema, rng, barcode, outprefix, outputbarcode) -> LongMoleculeRecipe:
    '''
    Randomly generates a long molecule and writes it to a FASTA file.
    Length of molecules is randomly distributed using an exponential distribution, with a minimum of 650bp.
    Returns a LongMoleculeRecipe that contains all the necessary information to simulate reads from that molecule.
    '''
    len_interval = schema.end - schema.start
    # make sure to cap the molecule length to the length of the interval/chromosome
    molecule_length = 0
    # make sure the molecule is greater than 650
    while molecule_length < 650:
        molecule_length = min(int(rng.exponential(scale = schema.mol_length)), len_interval)
    
    # set the max position to be length - mol_length to avoid additional computation
    adjusted_end = len_interval - molecule_length
    start = int(rng.uniform(low = 0, high = adjusted_end))
    end = start + molecule_length - 1

    molnumber = getrandbits(32)
    fasta_seq = schema.sequence[start-1:end+1]
    normalized_length = len(fasta_seq)-fasta_seq.count('N')

    if schema.mol_coverage < 1:
        N = max(1, normalized_length * schema.mol_coverage)/(schema.read_length*2)
    else:
        # draw N from a lognormal to avoid negative numbers, then convert the number back out of log space
        N = int(np.log(rng.lognormal(schema.mol_coverage, schema.mol_coverage/3)))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, normalized_length/(schema.read_length*2))
    if schema.singletons > 0:
        if rng.uniform(0,1) > schema.singletons:
            N = 1

    return LongMoleculeRecipe(schema.chrom, start, end, barcode, outputbarcode, molnumber, int(N), outprefix)
