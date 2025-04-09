#! /usr/bin/env python3

import os
import re
import sys
from datetime import datetime
import numpy as np
from rich.console import Console

mimick_console = Console(stderr=True, log_path=False)
RNG = np.random.default_rng()

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
	return [l[i:i+n] for i in range(0, len(l), n)]
