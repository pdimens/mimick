#! /usr/bin/env python3

import gzip
import multiprocessing
import os
import queue
import re
import shutil
import sys
from itertools import product
import pysam
from .common import error_terminate, mimick_console
from .classes import Schema, BarcodeGenerator
from .long_molecule import LongMoleculeRecipe

def readfq(fp):
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

def BEDtoInventory(bedfile, fasta, coverage, mol_cov, mol_len, read_len, singletons, circular) -> dict:
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
                inventory[idx] = Schema(haplotype, chrom,start, end, read_len, mean_reads_per, reads_required, mol_len, mol_cov, singletons, circular, _seq)
                idx += 1
    return inventory

def FASTAtoInventory(fasta, coverage, mol_cov, mol_len, read_len, singletons, circular) -> dict:
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
                    error_terminate(f"Error in {os.path.basename(fasta)} [yellow]contig {chrom}[/]: contigs must have at least 650 non-ambiguous ([yellow]N[/]) bases.")
                reads_required = int((coverage*normalized_length/read_len)/2)
                inventory[idx] = Schema(haplotype, chrom,start,end,read_len, mean_reads_per, reads_required, mol_len, mol_cov, singletons, circular, contig.sequence)
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
                error_terminate(f'Barcodes can only contain nucleotides [blue]A,T,C,G,U[/], but invalid barcode(s) provided: [yellow]{bc}[/]. This was first invalid barcode identified, but it may not be the only one.')

def interpret_barcodes(infile, lr_type) -> BarcodeGenerator:
    """
    Takes an open file connection and reads it line by line. Performs barcode validations and returns:
    - either an iter() or generator of barcodes (to use with next())
    - the total barcode	length (int)
    - the total number of barcodes [or combinations] (int)
    """
    try:
        with gzip.open(infile, 'rt') as filein:
            bc = list(set(i.strip() for i in filein.read().splitlines()))
    except gzip.BadGzipFile:
        with open(infile, 'r') as filein:
            bc = list(set(i.strip() for i in filein.read().splitlines()))
    except:
        error_terminate(f'Cannot open [yellow]{os.path.relpath(infile)}[/] for reading because it\'s not recognized as either a plaintext or gzipped file.')

    validate_barcodes(bc)
    bc_len = len(bc[0]) 
    if lr_type == "haplotagging":
        return BarcodeGenerator(product(bc,bc,bc,bc), lr_type, 2*bc_len, len(bc)**4)
    if lr_type == "stlfr":
        return BarcodeGenerator(product(bc,bc,bc,bc), lr_type, 3*bc_len, len(bc)**3)
    else:
        return BarcodeGenerator(iter(bc), lr_type, bc_len, len(bc))

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
        sequence = [f'@{name}#{outbc} {fr}', seq, '+', qual]
    return "\n".join(sequence)

class MoleculeRecorder:
    def __init__(self, prefix):
        self.inventory = open(f'{prefix}.molecules', 'w')
        self.inventory.write(
        "\t".join(["haplotype", "chromosome", "start_position", "end_position", "length", "reads", "nucleotide_barcode", "output_barcode"]) + "\n"
        )

    def close(self):
        self.inventory.close()

    def write(self, long_molecule: LongMoleculeRecipe):
        self.inventory.write(
            "\t".join([
                f"haplotype_{long_molecule.haplotype}",
                long_molecule.chrom,
                str(long_molecule.start),
                str(long_molecule.end),
                str(long_molecule.length),
                str(long_molecule.read_count),
                long_molecule.barcode,
                long_molecule.output_barcode
                ]) + "\n"
        )

class FileProcessor:
    def __init__(self, outprefix, outformat, quiet):
        self.output_format = outformat
        self.output_dir = os.path.dirname(outprefix)
        os.makedirs(self.output_dir, exist_ok=True)
        self.prefix = os.path.basename(outprefix)
        self.R1 = open(f"{outprefix}.R1.fq", 'w')
        self.R2 = open(f"{outprefix}.R2.fq", 'w')
        self.GFF = open(f"{outprefix}.gff", 'w')
        self.quiet = quiet == 2
        self.task_queue = multiprocessing.Queue()
        self.process = None
        self.running = False
        for ext in ["R1.fq", "R2.fq", "gff"]:
            if os.path.exists(f"{outprefix}.{ext}.gz"):
                os.remove(f"{outprefix}.{ext}.gz")

    def _worker(self):
        """Worker function that runs in the separate process"""
        while True:
            try:
                # Get task from queue (blocks until available)
                task = self.task_queue.get()
                # Shutdown signal
                if task is None:
                    # R1 read
                    if not self.quiet:
                        mimick_console.log(f"Compressing [blue]{os.path.basename(self.R1.name)}[/]")
                    pysam.tabix_compress(self.R1.name, f'{self.R1.name}.gz', force=True)
                    os.remove(self.R1.name)
                    # R2 read
                    if not self.quiet:
                        mimick_console.log(f"Compressing [blue]{os.path.basename(self.R2.name)}[/]")
                    pysam.tabix_compress(self.R2.name, f'{self.R2.name}.gz', force=True)
                    os.remove(self.R2.name)
                    # GFF file
                    if not self.quiet:
                        mimick_console.log(f"Compressing [blue]{os.path.basename(self.GFF.name)}[/]")
                    with open(self.GFF.name, 'rb') as f_in, gzip.open(f"{self.GFF.name}.gz", 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(self.GFF.name)
                    break

                _basename, barcode, output_barcode = task

                _basename = os.path.join(self.output_dir, "temp", _basename)
                try:
                    for i,f in enumerate([self.R1, self.R2],1):
                        with open(f"{_basename}.R{i}", 'r') as src:
                            for name,seq,qual in readfq(src):
                                if i== 2 and self.output_format == "10x":
                                    sequence = '\n'.join([f'@{name}',seq,'+',qual])
                                else:
                                    sequence = format_linkedread(
                                        name = name,
                                        bc = barcode,
                                        outbc = output_barcode,
                                        outformat = self.output_format,
                                        seq = seq,
                                        qual = qual,
                                        forward = i == 1
                                    )
                                    f.write(sequence + '\n')
                    os.remove(f"{_basename}.R{i}")

                    with open(f"{_basename}.gff", 'r') as src:
                        shutil.copyfileobj(src, self.GFF)
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
        self.task_queue.put((long_molecule.output_basename, long_molecule.barcode, long_molecule.output_barcode))

    def close_files(self):
        self.R1.close()
        self.R2.close()
        self.Gff.close()

    def error(self):
        """Stop the file processing cold"""
        try:
            self.close_files()
        except:
            pass
        if self.running:
            # Send shutdown signal
            self.task_queue.put(False)
            if self.process.is_alive():
                self.process.terminate()
            self.running = False

    def stop(self):
        """Stop the file processing process"""
        try:
            self.close_files()
        except:
            pass
        if self.running:
            # Send shutdown signal
            self.task_queue.put(None)
            # Wait for jobs to finish
            self.process.join()
            if self.process.is_alive():
                self.process.terminate()
            self.running = False
