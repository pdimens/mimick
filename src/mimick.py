#!/usr/bin/python3 env

import os
import sys
import glob
import re
import math
import gzip
import pysam
import multiprocessing
from itertools import product
from shutil import which

import pybedtools
from .common import *

def run(parser,args):

	'''
	Check arguments, run functions
	'''
	 # block pywgsim stdout
	redirect_stdout()
	print(f'[{get_now()}] mimick v{__version__}', file = sys.stderr)

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.BED=os.path.abspath(args.bedfile)
	c.FASTADIR=os.path.abspath(args.fasta)
	c.PREFIX=args.prefix
	c.threads=args.threads

	os.makedirs(c.OUT, exist_ok= True)
	if not os.access(os.path.abspath(c.OUT),os.W_OK):

		print(f'[{get_now()}][Error] Missing write permissions on the output folder', file = sys.stderr)
		sys.exit(1)

	if which('bedtools') is None:

		print(f'[{get_now()}][Error] bedtools must be in PATH', file = sys.stderr)
		sys.exit(1)

	try:

		bedfile=pybedtools.BedTool(c.BED)
		bedsrtd=bedfile.sort()

	except:

		print(f'[{get_now()}][Error] BED {c.BED} does not exist, is not readable, or is not a valid BED', file = sys.stderr)
		sys.exit(1)

	# read in the barcodes, then parse and validate
	c.barcodepath=args.barcodes
	c.barcodetype=args.barcode_type.lower()

	print(f'[{get_now()}] Validating the supplied barcodes', file = sys.stderr)
	try:
		with gzip.open(c.barcodepath, 'rt') as filein:
			c.barcodes, c.barcodebp, c.totalbarcodes = interpret_barcodes(filein, c.barcodetype)
	except gzip.BadGzipFile:
		with open(c.barcodepath, 'r') as filein:
			c.barcodes, c.barcodebp, c.totalbarcodes = interpret_barcodes(filein, c.barcodetype)
	except:
		print(f'[{get_now()}][Error] Cannot open {c.barcodepath} for reading', file = sys.stderr)
		sys.exit(1)

	# fill c with wgsim and general linked-read parameters 
	c.remainingbarcodes = c.totalbarcodes
	c.coverage=args.coverage
	c.error=args.error
	c.distance=args.distance
	c.stdev=args.stdev
	c.length=args.length
	c.mutation=args.mutation
	c.indels=args.indels
	c.extindels=args.extindels
	c.molnum=args.molecule_number
	c.mollen=args.molecule_length
	c.molcov=args.molecule_coverage
	c.barcodelist= open(f'{c.OUT}/{c.PREFIX}.barcodes', 'w')

	if c.barcodetype in ["10x", "tellseq"]:
		# barcode at beginning of read 1
		c.len_r1 = c.length - c.barcodebp
		if c.len_r1 <= 5:
			print(f'[{get_now()}][Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', file = sys.stderr)
			sys.exit(1)
		c.len_r2 = c.length
	elif c.barcodetype == "stlfr":
		# barcode at the end of read 2
		c.len_r1 = c.length
		c.len_r2 = c.length - c.barcodebp
		if c.len_r2 <= 5:
			print(f'[{get_now()}][Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', file = sys.stderr)
			sys.exit(1)
	else:
		# would be 4-segment haplotagging where AC on read 1 and BD on read 2
		c.len_r1 = c.length - c.barcodebp
		c.len_r2 = c.length - c.barcodebp
		if c.len_r1 <= 5 or len_r2 <= 5:
			print(f'[{get_now()}][Error] Removing barcodes from the reads would leave sequences <= 5 bp long. Read length: {c.length}, Barcode length: {c.barcodebp}', file = sys.stderr)
			sys.exit(1)
	if args.output_format:
		c.outformat=args.output_format.lower()
	else:
		c.outformat = c.barcodetype
	if c.outformat == "haplotagging":
		bc_range = [f"{i}".zfill(2) for i in range(1,97)]
		c.bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
	if c.outformat == "stlfr":
		bc_range = range(1, 1537)
		c.bc_generator = product(bc_range, bc_range, bc_range)

	fasta_files = [
		f for f in glob.glob(f'{os.path.abspath(c.FASTADIR)}/*') 
		if re.search(r'\.(fa|fasta)$', f, re.IGNORECASE)
	]
	if len(fasta_files) == 0:
		print(f'[{get_now()}][Error] No FASTA files detected in {c.FASTADIR}. If your FASTA files are gzipped, please decompress them.', file = sys.stderr)
		sys.exit(1)
	c.ffiles=sorted(fasta_files, key=natural_keys) #list all FASTA in folder
	c.regioncoverage=c.coverage/len(c.ffiles)
	
	# check that the haplotagging output format can support the supplied number of barcodes
	if c.outformat == "haplotagging":
		if c.totalbarcodes > 96**4:
			print(f'[{get_now()}][Error] The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in haplotagging format is limited to {96**4} barcodes', file = sys.stderr)
			sys.exit(1)
	# check that the stlfr output format can support the supplied number of barcodes
	if c.outformat == "stlfr":
		if c.totalbarcodes > 1537**3:
			print(f'[{get_now()}][Error] The barcodes and barcode type supplied will generate a potenial {c.totalbarcodes} barcodes, but outputting in stlfr format is limited to {1537**3} barcodes', file = sys.stderr)
			sys.exit(1)

	print(f'[{get_now()}] Preparing for bulk simulations with a single clone', file = sys.stderr) 

	for k,s in enumerate(c.ffiles):

		print(f'[{get_now()}] Processing haplotype {k+1}', file = sys.stderr)
		c.hapnumber=f'{k+1}'
		c.ffile=c.ffiles[k]

		for w in bedsrtd:
	
			print(f'[{get_now()}] Simulating from region {w.chrom}:{w.start}-{w.end}', file = sys.stderr)
			LinkedSim(w,c)

	allfastq=glob.glob(os.path.abspath(c.OUT) + '/*.fastq')

	#gzip multiprocessing

	chunk_size=len(allfastq)/c.threads
	slices=Chunks(allfastq,math.ceil(chunk_size))
	processes=[]

	for i,sli in enumerate(slices):
		for _s in sli:
			print(f'[{get_now()}] Compressing {os.path.basename(_s)}', file = sys.stderr)

		p=multiprocessing.Process(target=BGzipper, args=(sli,))
		p.start()
		processes.append(p)
		
	for p in processes:
		
		p.join()
	c.barcodelist.close()
	print(f'[{get_now()}] Done', file = sys.stderr)
