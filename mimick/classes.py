#! /usr/bin/env python3

import os
from itertools import product

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

    def adjust_readlength(self, lr_chem, bc_bp):
        """
        Modifiest the lenth of the R1 and R2 read lengths based on the linked-read chemistry
        """
        if lr_chem in ["10x", "tellseq"]:
            # barcode at beginning of read 1
            if self.length_R1 - bc_bp <= 5:
                error_terminate(f'Removing barcodes from the reads would leave R1 sequences <= 5 bp long. Read length: {self.length_R1}, Barcode length: {bc_bp}')
            self.length_R1 -= bc_bp
        elif lr_chem == "stlfr":
            # barcode at the end of read 2
            if self.length_R2 - bc_bp <= 5:
                error_terminate(f'Removing barcodes from the reads would leave R2 sequences <= 5 bp long. Read length: {self.length_R1}, Barcode length: {bc_bp}')
            self.length_R2 -= bc_bp

class BarcodeGenerator():
    '''
    The container for generating barcodes on the fly, where bc_in is the original nucleotide barcodes and bc_out are the output (translated) barcode
    '''
    def __init__(self, generator, lr_type, bc_len, bc_total):
        self.lr_type = lr_type
        self.bc_in = generator
        self.max = bc_total
        self.bc_length = bc_len
        self.remaining = 0
    
    def setup_output_bc(self, lr_type):
        """
        Based on the output linked-read type, sets up a barcode generator. If 10X/tellseq, doesn't set up a generator b/c one isn't needed.
        """
        if lr_type == "haplotagging":
            if self.max > 96**4:
                error_terminate(f'The barcodes and barcode type supplied will generate a potential {self.max} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')
            bc_range = [f"{i}".zfill(2) for i in range(1,97)]
            self.bc_out_type = "haplotagging"
            self.bc_out = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
        elif lr_type == "stlfr":
            if self.max > 1537**3:
                error_terminate(f'The barcodes and barcode type supplied will generate a potential {self.max} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes')
            bc_range = [str(i) for i in range(1, 1537)]
            self.bc_out_type = "stlfr"
            self.bc_out = product(bc_range, bc_range, bc_range)
        else:
            self.bc_out_type = "tellseq"
            self.bc_out = None

    def get_next_out(self, nucleotides = None):
        """
        If self.lr_type is haplotagging or stlfr, returns a barcode of that style,
        otherwise returns the nucleotide barcode used as input
        """
        if self.bc_out_type == "haplotagging":
            _char = ""
        elif self.bc_out_type == "stlfr":
            _char = "_"
        else:
            return nucleotides
        return _char.join(next(self.bc_out))
    
    def get_next_bc(self):
        """
        Properly format the next nucleotide barcode
        """
        return "".join(next(self.bc_in))

    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            outstring += f"{i}: {j}\n"
        return outstring

class Schema():
    '''
    The container for per-chromosome/contig parameters. This is what guides
    the creation of LongMolecules. Values like read_length/read_pairs_per_molecule,
    mol_length etc. are averages. A contig/interval being circular will have its
    sequence repeated once, but the start and end positions don't change
    '''
    def __init__(self, haplotype, chrom, start, end, read_length, read_pairs_per_mol, reads_required, mol_length, mol_cov, singletons, circular, seq):
        self.haplotype = haplotype
        self.chrom = chrom
        self.start = start
        self.end = end
        self.read_length = int(read_length)
        self.read_pairs_per_mol = int(read_pairs_per_mol)
        self.reads_current = 0
        self.reads_required = reads_required
        self.mol_length = mol_length
        self.mol_coverage = mol_cov
        self.singletons = singletons
        self.is_circular = circular
        if circular:
            self.sequence = seq * 2
        else:
            self.sequence = seq

    def __str__(self):
        outstring = ""
        for i,j in self.__dict__.items():
            if i != "sequence":
                outstring += f"{i}: {j}\n"
            else:
                outstring += f"{i}: {j[:31]}... (length = {self.end+1 - self.start})"
        return outstring
