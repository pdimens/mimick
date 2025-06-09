#! /usr/bin/env python3

import os
import re
import sys
from itertools import product
import pysam
from .common import error_terminate, mimick_console
from .classes import Schema

def readfq(fp): # this is a fast generator function
    '''Yield FASTQ record'''
    read=[]
    for line in fp:
        read.append(line.rstrip())
        if len(read) == 4:
            yield read[0][1:], read[1], read[3]
            read=[]

def index_fasta(fasta):
    try:
        _fa = os.path.abspath(fasta)
        pysam.faidx(_fa)
        if _fa.lower().endswith(".gz") and not os.path.exists(_fa + ".gzi"):
            os.system(f"bgzip --reindex {_fa}")
    except Exception as e:
        error_terminate(f'Failed to index {_fa}. Error reported by samtools:\n{e}', False)

def BEDtoInventory(bedfile, fasta, coverage, mol_cov, mol_len, read_len, singletons) -> dict:
    '''
    Read the BED file, do validation against the FASTA files and derive the schema, and return a dict of Schema objects that's an
    inventory tracker in the form of d[idx] = [read_count, reads_requireduired, Schema]
    '''
    inventory = {}
    mean_reads_per = (mol_cov*mol_len)/(read_len*2) if mol_cov < 1 else mol_cov
    mean_reads_per = max(1, mean_reads_per)
    haplotype = 0
    idx = 0
    for _fasta in fasta:
        haplotype += 1
        with open(bedfile, "r", encoding="utf-8") as bed, pysam.FastxFile(_fasta) as fa:
            for line in bed:
                row = line.split()
                try:
                    chrom = row[0]
                    start = int(row[1])
                    end = int(row[2])
                except ValueError:
                    error_terminate(f"The input file is formatted incorrectly at line {idx+1}. This is the first row triggering this error, but it may not be the only one.")
                if start > end:
                    error_terminate(f"The interval start position is greater than the interval end position at line {idx+1}. This is the first row triggering this error, but it may not be the only one.")
                if (end - start) < 650:
                    error_terminate(f"Error in {os.path.basename(bedfile)} [yellow]line {idx+1}[/]: the designated interval must be at least 650bp long.")

                _seq = _fasta.fetch(chrom, start-1, end+1)
                normalized_length = (end-start) - _seq.count('N')
                reads_required = int((coverage*normalized_length/read_len)/2)
                inventory[idx] = Schema(haplotype, chrom,start, end, read_len, mean_reads_per, reads_required, mol_len, mol_cov, singletons, _seq)
                idx += 1
    return inventory

def FASTAtoInventory(fasta, coverage, mol_cov, mol_len, read_len, singletons) -> dict:
    '''
    Read the FASTA files and derive the contig name, start, and end positions and other simulation schema
    and return a dict of Schema objects that's an inventory tracker in the form of
    d[idx] = Schema
    '''
    inventory = {}
    mean_reads_per = (mol_cov*mol_len)/(read_len*2) if mol_cov < 1 else mol_cov
    mean_reads_per = max(1, mean_reads_per)
    haplotype = 0
    idx = 0
    for _fasta in fasta:
        haplotype += 1
        with pysam.FastxFile(_fasta) as fa:
            for contig in fa:
                chrom = contig.name
                start = 1
                end = len(contig.sequence)
                normalized_length = end - contig.sequence.count('N')
                if normalized_length < 650:
                    error_terminate(f"Error in {os.path.basename(fasta)} [yellow]contig {chrom}[/]: contigs must have at least 650 non-ambiguous (N) bases.")
                reads_required = int((coverage*normalized_length/read_len)/2)
                inventory[idx] = Schema(haplotype, chrom,start,end,read_len, mean_reads_per, reads_required, mol_len, mol_cov, singletons, contig.sequence)
                idx += 1
    return inventory

def validate_barcodes(bc_list):
    '''Takes a file with barcodes and validates them to be ATGCU nucleotides and barcodes same length'''
    # check first row for multiple columns, if there are multiple, it's haplotagging
    if len(bc_list[0].strip().split()) != 1:
        error_terminate(f'Barcode file is expected to only have one barcode per line')
    else:
        bc_lens = set()
        for i in bc_list:
            bc_lens.add(len(i))
            if len(bc_lens) > 1:
                error_terminate(f'Barcodes provided must all be the same length')
        # validate barcodes are only ATCGU nucleotides
        for bc in bc_list:
            if not bool(re.fullmatch(r'^[ATCGU]+$', bc, flags = re.IGNORECASE)):
                error_terminate(f'Barcodes can only contain nucleotides A,T,C,G,U, but invalid barcode(s) provided: {bc}. This was first invalid barcode identified, but it may not be the only one.')
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
