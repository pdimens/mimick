#! /usr/bin/env python3

import multiprocessing
import queue
import shutil
import time
from pathlib import Path
import os
import re
import sys
from itertools import product
import pysam
from .common import error_terminate, mimick_console
from .classes import Schema
from .long_molecule import LongMoleculeRecipe

def readfq(fp): # this is a fast generator function
    '''Yield FASTQ record'''
    read=[]
    for line in fp:
        read.append(line.rstrip())
        if len(read) == 4:
            yield read[0][1:], read[1], read[3]
            read=[]

def index_fasta(fasta):
    outs = []
    try:
        _fa = os.path.abspath(fasta)
        pysam.faidx(_fa)
        outs.append(_fa + ".fai")
        if _fa.lower().endswith(".gz") and not os.path.exists(_fa + ".gzi"):
            os.system(f"bgzip --reindex {_fa}")
            outs.append(_fa + ".gzi")
    except Exception as e:
        error_terminate(f'Failed to index {_fa}. Error reported by samtools:\n{e}', False)
    return outs

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

def format_linkedread(name, bc, outbc, outformat, seq, qual, forward: bool):
    '''Given a linked-read output type, will format the read accordingly and return it'''
    if outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * len(bc)}{qual}']
    elif outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name}:{bc} {fr}', seq, '+', qual]
    elif outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        sequence = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{outbc}', seq, '+', qual]
    elif outformat == "standard":
        fr = "/1" if forward else "/2"
        sequence = [f'@{name}{fr}\tVX:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    else:
        print(outformat)
        sys.exit(1)
    return "\n".join(sequence)

class FileProcessor:
    def __init__(self, outprefix, outformat, quiet):
        self.output_format = outformat
        self.output_dir = os.path.dirname(outprefix)
        os.makedirs(self.output_dir, exist_ok=True)
        self.prefix = os.path.basename(outprefix)
        self.R1 = f"{outprefix}.R1.fq"
        self.R2 = f"{outprefix}.R2.fq"
        self.GFF = f"{outprefix}.gff"
        self.quiet = quiet == 2
        self.task_queue = multiprocessing.Queue()
        self.process = None
        self.running = False

    def _worker(self):
        """Worker function that runs in the separate process"""
        while True:
            try:
                # Get task from queue (blocks until available)
                #task = self.task_queue.get(timeout=1)
                task = self.task_queue.get()
                if task is None:  # Shutdown signal
                    if not self.quiet:
                        mimick_console.log(f"Compressing [blue]{os.path.basename(self.R1)}[/]")
                    pysam.tabix_compress(self.R1, f'{self.R1}.gz', force=True)
                    os.remove(self.R1)
                    if not self.quiet:
                        mimick_console.log(f"Compressing [blue]{os.path.basename(self.R2)}[/]")
                    pysam.tabix_compress(self.R1, f'{self.R2}.gz', force=True)
                    os.remove(self.R2)
                    break
                if task is False:
                    # Exit on error, don't compress
                    self.task_queue.task_done()
                    break
                _basename, barcode, output_barcode = task

                _basename = os.path.join(self.output_dir, "temp", _basename)
                try:
                    with open(self.R1, 'a') as out, open(f"{_basename}.R1", 'r') as src:
                        for name,seq,qual in readfq(src):
                            sequence = format_linkedread(
                                name = name,
                                bc = barcode,
                                outbc = output_barcode,
                                outformat = self.output_format,
                                seq = seq,
                                qual = qual,
                                forward = True
                            )
                            out.write(sequence + '\n')
                    
                    os.remove(f"{_basename}.R1")

                    with open(self.R2, 'a') as out, open(f"{_basename}.R2", 'r') as src:
                        for name,seq,qual in readfq(src):
                            if self.output_format == "10x":
                                sequence = '\n'.join([f'@{name}',seq,'+',qual])
                            else:
                                sequence = format_linkedread(
                                    name = name,
                                    bc = barcode,
                                    outbc = output_barcode,
                                    outformat = self.output_format,
                                    seq = seq,
                                    qual = qual,
                                    forward = False
                                )
                            out.write(sequence + '\n')

                    os.remove(f"{_basename}.R2")

                    with open(self.GFF, 'a') as dest, open(f"{_basename}.gff", 'r') as src:
                        shutil.copyfileobj(src, dest)
                    
                    os.remove(f"{_basename}.gff")
                    
                except Exception as e:
                    self.process.terminate()
                    error_terminate(f"{e}")
                    
            except queue.Empty:
                continue
            except Exception as e:
                self.process.terminate()
                error_terminate(f"{e}")

    def start(self):
        """Start the file processing process"""
        if not self.running:
            self.process = multiprocessing.Process(target=self._worker)
            self.process.start()
            self.running = True

    def submit_files(self, long_molecule: LongMoleculeRecipe):
        """Submit temp files for processing"""
        if not self.running:
            raise RuntimeError("FileProcessor not started. Call start() first.")
        _basename = long_molecule.output_basename
        barcode = long_molecule.barcode
        output_barcode = long_molecule.output_barcode
        self.task_queue.put((_basename, barcode, output_barcode))

    def stop(self):
        """Stop the file processing process"""
        if self.running:
            self.task_queue.put(None)  # Send shutdown signal
            self.process.join()
            if self.process.is_alive():
                self.process.terminate()
            self.running = False
