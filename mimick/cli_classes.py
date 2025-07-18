#! /usr/bin/env python3

import click as _click

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