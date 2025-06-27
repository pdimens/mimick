#! /usr/bin/env python3

import os
import click as _click

class Barcodes(_click.ParamType):
    """A class for a click type which accepts either a file or two integers, separated by a comma."""
    name = "barcodes"
    def convert(self, value, param, ctx):
        if os.path.isfile(value):
            return os.path.abspath(value)
        try:
            bp,count = value.split(",")
        except ValueError:
            self.fail(f"{value} is not a file, nor in int,int format", param, ctx)
        try:
            bp = int(bp)
        except ValueError:
            self.fail(f"{bp} is not an integer.", param, ctx)
        if bp <= 2:
            self.fail(f"Barcodes must be at least 3 base pairs long. (input: {bp})", param, ctx)
        try:
            count = int(count)
        except ValueError:
            self.fail(f"{count} is not an integer.", param, ctx)
        if count <= 2:
            self.fail(f"There must be at least 3 barcodes. (input: {count})", param, ctx)
        return [bp, count]

class ReadLengths(_click.ParamType):
    """A class for a click type which accepts two integers, separated by a comma."""
    name = "readlengths"
    def convert(self, value, param, ctx):
        try:
            R1,R2 = value.split(",")
        except ValueError:
            self.fail(f"{value} is not in int,int format", param, ctx)
        try:
            R1 = int(R1)
        except ValueError:
            self.fail(f"{R1} is not an integer.", param, ctx)
        try:
            R2 = int(R2)
        except ValueError:
            self.fail(f"{R2} is not an integer.", param, ctx)
        if R1 <10:
            self.fail(f"R1 reads must must be at least 10 bp. (input: {R1})", param, ctx)
        if R2 <10:
            self.fail(f"R2 reads must be at least 10 bp. (input: {R2})", param, ctx)
        return [R1, R2]