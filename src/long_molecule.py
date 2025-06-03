#! /usr/bin/env python3

import os
import re
from random import getrandbits
from .classes import Schema

class LongMoleculeRecipe(object):
    '''Molecule instance'''
    def __init__(self,fasta,barcode,outbarcode,chrom,read_count,mol_id):
        self.fasta=fasta
        self.barcode=barcode
        self.output_barcode=outbarcode
        self.chrom=chrom
        self.read_count=read_count
        self.mol_id=mol_id
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            if i != "seq":
                outstring += f"{i}: {j}\n"
            else:
                outstring += f"{i}: " + j[:min(30, len(j))] + f"(length = {len(j)})\n"
        return outstring

def create_long_molecule(schema: Schema, rng, barcode, outprefix, outputbarcode) -> LongMoleculeRecipe:
    '''
    Randomly generates a long molecule and writes it to a FASTA file.
    Length of molecules is randomly distributed using an exponential distribution, with a minimum of 300bp.
    Returns a LongMoleculeRecipe that contains all the necessary information to simulate reads from that molecule.
    '''
    lensingle = len(schema.seq)
    start = int(rng.uniform(low = 0, high = lensingle))
    while True:
        # keep trying until a molecule of at least 300bp is created
        length = int(rng.exponential(scale = schema.mol_length))
        end = min(start + length - 1, lensingle)
        if end - start >= 300:
            break
    molnumber = getrandbits(32)
    fasta_header = f'>CHROM:{schema.chrom}_START:{start}_END:{end}_BARCODE:{barcode}_MOL:{molnumber}'
    fasta_seq = schema.seq[start-1:end+1]
    molecule_fasta = os.path.abspath(f'{outprefix}{schema.chrom}.{barcode}.{molnumber}.fa')
    with open(molecule_fasta, 'w') as faout:
        faout.write(
            "\n".join([fasta_header, '\n'.join(re.findall('.{1,60}', fasta_seq))]) + "\n"
        )

    normalized_length = len(fasta_seq)-fasta_seq.count('N')
    if schema.mol_coverage < 1:
        N = int(normalized_length * schema.mol_coverage)/(schema.read_length*2)
    else:
        # draw N from a normal distribution with a mean of molcov and stdev of molcov/3, avoiding < 0
        N = max(0, rng.normal(schema.mol_coverage, schema.mol_coverage/3))
        # set ceiling to avoid N being greater than can be sampled
        N = min(N, normalized_length/(schema.read_length*2))
    if schema.singletons > 0 and N != 0:
        if rng.uniform(0,1) > schema.singletons:
            N = 1    
    return LongMoleculeRecipe(molecule_fasta, barcode, outputbarcode, schema.chrom, int(N), molnumber)