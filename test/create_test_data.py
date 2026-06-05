#! /usr/bin/env python3
"""Creates a FASTA file with 2 contigs of 200k random nucleotides and 1 contig of 1k nucleotides"""
from random import choices
import gzip
import os
import sys

prefix = sys.argv[1] if len(sys.argv) > 1 else "."

def fasta():
    nuc = "ATCG"
    return (
        f">Contig1\n{''.join(choices(nuc, k = 200000))}\n"
        f">Contig2\n{''.join(choices(nuc, k = 200000))}\n"
        f">Contig3\n{''.join(choices(nuc, k = 1000))}NNNN\n"
    ).encode()

for i in [1,2]:
    with gzip.open(os.path.join(prefix, f"hap{i}.fa.gz"), "wb") as fa:
        _ = fa.write(fasta())
