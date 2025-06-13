#! /usr/bin/env python3

import os
from pywgsim import wgsim
from .classes import wgsimParams
from .common import *
from .long_molecule import LongMoleculeRecipe

def linked_simulation(wgsimparams: wgsimParams, long_molecule: LongMoleculeRecipe):
    '''
    Process a LongMoleculeRecipe into reads using [py]wgsim. Returns the input LongMoleculeRecipe to be parsed
    later in the quota/progressbar and final output appending process.
    '''
    tempdir = os.path.join(wgsimparams.outdir,"temp")
    R1 = f"{tempdir}/{long_molecule.output_basename}.R1"
    R2 = f"{tempdir}/{long_molecule.output_basename}.R2"
    gff = f"{tempdir}/{long_molecule.output_basename}.gff"
    try:
        wgsim.core(
            long_molecule.fasta,
            r1 = R1,
            r2 = R2,
            err_rate = wgsimparams.error,
            mut_rate = wgsimparams.mutation,
            indel_frac = wgsimparams.indels,
            indel_ext = wgsimparams.extindels,
            max_n = 0.05,
            is_hap = 0,
            N = long_molecule.read_count,
            dist = wgsimparams.read_distance,
            stdev = wgsimparams.distance_stdev,
            size_l = wgsimparams.length_R1,
            size_r = wgsimparams.length_R2,
            seed = wgsimparams.randomseed,
            gff = gff
        )

    except KeyboardInterrupt:
        mimick_keyboardterminate()
    except Exception as e:
        return f"{e}"
    finally:
        os.remove(long_molecule.fasta)

    if os.stat(R1).st_size == 0:
        os.remove(R1)
        os.remove(R2)
        os.remove(gff)
        return 0

    return long_molecule
