#! /usr/bin/env python3

import os

class wgsimParams():
    '''
    The container for pywgsim parameters. This gets sent per-thread for execution.
    '''
    def __init__(self, error, mutation, indels, extindels, distance, stdev, length_r1, length_r2, seed, outprefix):
        self.error = error
        self.mutation = mutation
        self.indels = indels
        self.extindels = extindels
        self.read_distance = distance
        self.distance_stdev = stdev
        self.length_R1 = length_r1
        self.length_R2 = length_r2
        self.randomseed = seed
        self.outdir =  os.path.dirname(outprefix)
        self.prefix = os.path.basename(outprefix)
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            outstring += f"{i}: {j}\n"
        return outstring

class Schema():
    '''
    The container for per-chromosome/contig parameters. This is what guides
    the creation of LongMolecules. Values like read_length/read_pairs_per_molecule,
    mol_length etc. are averages.
    '''
    def __init__(self, haplotype, chrom, start, end, read_length, read_pairs_per_mol, reads_required, mol_length, mol_cov, singletons, seq):
        self.haplotype = haplotype
        self.chrom = chrom
        self.start = start
        self.end = end
        self.read_length = read_length
        self.read_pairs_per_mol = read_pairs_per_mol
        self.reads_current = 0
        self.reads_required = reads_required
        self.mol_length = mol_length
        self.mol_coverage = mol_cov
        self.singletons = singletons
        self.sequence = seq
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            if i != "sequence":
                outstring += f"{i}: {j}\n"
            else:
                outstring += f"{i}: {j[:31]}... (length = {len(j)})"
        return outstring
