#! /usr/bin/env python3

import os
import re
import sys
from datetime import datetime
from itertools import product
from random import sample
from rich.console import Console
from rich.table import Table
from rich import box

mimick_console = Console(stderr=True, log_path=False)

def mimick_keyboardterminate():
    mimick_console.print("")
    mimick_console.rule("[bold]Terminating Mimick", style = "yellow")
    sys.exit(1)

def error_terminate(text: str, rule: bool = True):
    mimick_console.log(f"[Error] {text}", highlight=False, style = "red")
    if rule:
        mimick_console.rule("[bold]Terminating Mimick due to an error", style = "red")
    sys.exit(1)

def log_table():
    '''Create a pre-formatted Rich Table to add rows to later'''
    _t = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), min_width = 55, box=box.SIMPLE)
    _t.add_column("detail", justify="left", no_wrap=True)
    _t.add_column("value", style="magenta", justify="right")
    return _t

def redirect_stdout():
    '''Suppress c/c++/fortran stdout, keep python print calls'''
    # flush is important when redirecting to files
    sys.stdout.flush()
    newstdout = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(newstdout, 'w')

def atoi(text):
    '''Convert text to integers'''
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''Natural sort'''
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def Chunks(l,n):
    '''Split list in chunks based on number of threads'''
    n = max(1,n)
    return [l[i:i+n] for i in range(0, len(l), n)]

def generate_barcodes(bp):
    '''
    Given a `bp` length, create a barcode generator. The nucleotide orders are randomized for each position
    so the barcodes don't all start with AAAAAAAAAAAAA and look more visually distinct. The randomization doesn't
    serve a functional purpose beyond that.
    '''
    return product(*[sample("ATCG", 4) for i in range(bp)])

#TODO this needs a significant rework
def format_linkedread(name, bc, outbc, outformat, seq, qual, forward: bool):
    '''Given a linked-read output type, will format the read accordingly and return it'''
    if outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * Container.barcodebp}{qual}']
    elif outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}:{bc} {fr}', seq, '+', qual]
    elif outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        read = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{outbc}', seq, '+', qual]
    elif outformat == "standard":
        fr = "/1" if forward else "/2"
        read = [f'@{name}{fr}\tVX:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    return read