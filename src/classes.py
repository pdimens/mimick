#! /usr/bin/env python3

class c():
    '''Container. This stores argparser parameters. Used to pass multiple parameters at once.'''
    OUT = ''
    PREFIX = ''
    FASTADIR = ''

    #pywgsim
    coverage=0
    regioncoverage=0
    error=0
    distance=0
    stdev=0
    length=0
    length_r1=0
    length_r2=0
    mutation=0
    indels=0
    extindels=0

    #bulk
    ffiles=None
    ffile=None
    outformat=None
    hapnumber=0
    threads=0

    #molecules
    barcodepath=None
    molnum=0
    mollen=0
    molcov=0
    barcodetype=None
    barcodebp=0
    barcodes=None	# will be an iterable to use with next()
    bc_generator = None
    used_bc={}
    totalbarcodes=0
    remainingbarcodes=0
    barcodeslist=None

class Interval():
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

class Molecule(object):
    '''Molecule instance'''
    def __init__(self,length,start,end,index):
        self.seqidx=index
        self.length=length
        self.index_droplet=0
        self.barcode=None
        self.start=start
        self.end=end
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)
