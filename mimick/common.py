#! /usr/bin/env python3

import sys
from time import gmtime, strftime
from itertools import product
from random import sample
from rich.console import Console
from rich.progress import Progress, TextColumn, TimeElapsedColumn, TaskProgressColumn, BarColumn

mimick_console = Console(stderr=True, log_path=False)
STYLES = ["purple", "yellow", "green", "orange", "blue", "magenta"] * 4

PROGRESS = Progress(
    TextColumn("[progress.description]{task.description}"),
    BarColumn(finished_style="purple", complete_style="yellow"),
    TaskProgressColumn(),
    TimeElapsedColumn(),
    transient=True,
    console=mimick_console
)

def mimick_keyboardterminate():
    '''Print text saying mimick is terminating and kill it'''
    mimick_console.print("")
    mimick_console.rule("[bold]Terminating Mimick", style = "yellow")
    PROGRESS.stop()
    sys.exit(1)

def error_terminate(text: str, rule: bool = True, appender = None):
    if rule:
        mimick_console.rule(f"Error", style = "red")
    fmt_time = strftime("Date: %Y-%m-%d Time: %H:%M:%S", gmtime())
    mimick_console.print(f"{text}", highlight=False, style = "red")
    if rule:
        mimick_console.rule(f"[dim]{fmt_time}[/]", style = "dim")
    PROGRESS.stop()
    if appender:
        appender.error()
    sys.exit(1)

def generate_barcodes(bp):
    '''
    Given a `bp` length, create a barcode generator. The nucleotide orders are randomized for each position
    so the barcodes don't all start with AAAAAAAAAAAAA and look more visually distinct. The randomization doesn't
    serve a functional purpose beyond that.
    '''
    return product(*[sample("ATCG", 4) for i in range(bp)])
