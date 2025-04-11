#! /usr/bin/env python3

import os
import re
import sys
from itertools import product
import click as _click
import pysam
from .common import mimick_console
from .classes import Interval

def readfq(fp): # this is a fast generator function
    '''Yield FASTQ record'''
    read=[]
    for line in fp:
        read.append(line.rstrip())
        if len(read) == 4:
            yield read[0][1:], read[1], read[3]
            read=[]

def BGzipper(sli,):
    '''Use pysam/htslib BGzip and multi-processing to save some time'''
    for s in sli:
        mimick_console.log(f'Compressing [blue]{os.path.basename(s)}[/]')
        pysam.tabix_compress(s, f'{s}.gz', force=True)
        os.remove(s)

def FASTAtoBED(fasta):
    '''Read a FASTA file and derive the contig name, start, and end positions'''
    intervals = []
    with pysam.FastxFile(fasta) as fa:
        for contig in fa:
            chrom = contig.name
            start = 1
            end = len(contig.sequence)
            intervals.append(Interval(chrom,start,end))
    return intervals

def readBED(bedfile):
    '''Read the BED file, do basic validation, and return a list of Interval objects'''
    intervals = []
    with open(bedfile, "r", encoding="utf-8") as bed:
        for idx, line in enumerate(bed, 1):
            row = line.split()
            try:
                start = int(row[1])
                end = int(row[2])
            except ValueError:
                mimick_console.log(f"[Error] The input file is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.", highlight=False, style = "red")
                sys.exit(1)
            if start > end:
                mimick_console.log(f"[Error] The interval start position is greater than the interval end position at row {idx}. This is the first row triggering this error, but it may not be the only one.", highlight=False, style = "red")
                sys.exit(1)
            intervals.append([row[0], start, end])
    return [Interval(*i) for i in sorted(intervals)]

def validate_barcodes(bc_list):
    '''Takes a file with barcodes and validates them to be ATGCU nucleotides and barcodes same length'''
    # check first row for multiple columns, if there are multiple, it's haplotagging
    if len(bc_list[0].strip().split()) != 1:
        mimick_console.log(f'[Error] Barcode file is expected to only have one barcode per line', highlight=False, style = "red")
        sys.exit(1)
    else:
        bc_lens = set()
        for i in bc_list:
            bc_lens.add(len(i))
            if len(bc_lens) > 1:
                mimick_console.log(f'[Error] Barcodes provided must all be the same length', highlight=False, style = "red")
                sys.exit(1)
        # validate barcodes are only ATCGU nucleotides
        for bc in bc_list:
            if not bool(re.fullmatch(r'^[ATCGU]+$', bc, flags = re.IGNORECASE)):
                mimick_console.log(f'[Error] Barcodes can only contain nucleotides A,T,C,G,U, but invalid barcode(s) provided: {bc}. This was first invalid barcode identified, but it may not be the only one.', highlight=False, style = "red")
                sys.exit(1)

def interpret_barcodes(infile, lr_type):
    """
    Takes an open file connection and reads it line by line. Performs barcode validations and returns:
    - either an iter() or generator of barcodes (to use with next())
    - the total barcode	length (int)
    - the total number of barcodes [or combinations] (int)
    """
    bc = list(set(i.strip() for i in infile.read().splitlines()))
    validate_barcodes(bc)
    bc_len = len(bc[0]) 
    if lr_type == "haplotagging":
        # 2 barcodes per
        return product(bc,bc,bc,bc), 2 * bc_len, len(bc)**4
    if lr_type == "stlfr":
        return product(bc,bc,bc), 3 * bc_len, len(bc)**3
    else:
        return iter(bc), bc_len, len(bc)

class Barcodes(_click.ParamType):
    """A class for a click type which accepts either a file or two integers, separated by a comma."""
    name = "barcodes"
    def convert(self, value, param, ctx):
        if os.path.isfile(value):
            return os.path.abspath(value)
        try:
            bp,count = value.split(",")
        except ValueError:
            self.fail(f"{value} is not a file, not in int,int format", param, ctx)
        try:
            bp = int(bp)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        try:
            count = int(count)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        return [bp, count]
