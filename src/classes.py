#! /usr/bin/env python3

import os

class wgsimParams():
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
    def __init__(self, chrom, start, end, length, read_pairs_per_mol, reads_req, n_mol, mol_length, mol_cov, seq):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.read_length = length
        self.read_pairs_per_mol = read_pairs_per_mol
        self.reads_req = reads_req
        self.n_mol = n_mol
        self.mol_length = mol_length
        self.seq = seq
        self.singletons = None
        self.mol_coverage = mol_cov
    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            if i != "seq":
                outstring += f"{i}: {j}\n"
            else:
                outstring += f"{i}: " + j[:min(30, len(j))] + f"(length = {len(j)})\n"
        return outstring

