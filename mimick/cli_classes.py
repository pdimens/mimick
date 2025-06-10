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
            self.fail(f"{value} is not a file, not in int,int format", param, ctx)
        try:
            bp = int(bp)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        try:
            count = int(count)
        except ValueError:
            self.fail(f"{value} is not an integer.", param, ctx)
        return [bp, count]
