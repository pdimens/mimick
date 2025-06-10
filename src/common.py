#! /usr/bin/env python3

import sys
from itertools import product
from random import sample
from rich.console import Console
from rich.progress import Progress, TextColumn, TimeElapsedColumn, TaskProgressColumn, BarColumn

mimick_console = Console(stderr=True, log_path=False)

PROGRESS = Progress(
    TextColumn("[progress.description]{task.description}"),
    BarColumn(finished_style="purple", complete_style="yellow"),
    TaskProgressColumn(),
    TimeElapsedColumn(),
    transient=True,
    console=mimick_console
)

def mimick_keyboardterminate():
    mimick_console.print("")
    mimick_console.rule("[bold]Terminating Mimick", style = "yellow")
    sys.exit(1)

def error_terminate(text: str, rule: bool = True):
    mimick_console.log(f"[Error] {text}", highlight=False, style = "red")
    if rule:
        mimick_console.rule("[bold]Terminating Mimick due to an error", style = "red")
    sys.exit(1)

def generate_barcodes(bp):
    '''
    Given a `bp` length, create a barcode generator. The nucleotide orders are randomized for each position
    so the barcodes don't all start with AAAAAAAAAAAAA and look more visually distinct. The randomization doesn't
    serve a functional purpose beyond that.
    '''
    return product(*[sample("ATCG", 4) for i in range(bp)])
