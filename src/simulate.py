#! /usr/bin/env python3

import os
import sys
import re
import math
import multiprocessing
import time
from tempfile import mkstemp
from .classes import *
from .file_ops import *
from .common import *
import pyfaidx
from pywgsim import wgsim

def wrap_wgsim(*args, **kwargs):
    '''prevents wgsim.core output from being printed to the terminal'''
    # Create temp files for stdout and stderr
    stdout_fd, stdout_path = mkstemp()
    # Save original file descriptors
    orig_stdout = os.dup(1)
    orig_stderr = os.dup(2)
    try:
        # Redirect stdout/stderr at OS level
        os.dup2(stdout_fd, 1)
        os.dup2(stdout_fd, 2)
        # Call the function
        wgsim.core(*args, **kwargs)
    finally:
        # Restore original stdout/stderr
        os.dup2(orig_stdout, 1)
        os.dup2(orig_stderr, 2)
        # Close the temp files
        os.close(stdout_fd)
        os.close(orig_stdout)
        os.close(orig_stderr)
    # Read the captured output
    with open(stdout_path, 'r') as f:
        stdout_content = f.read()
    # Clean up temp files
    os.unlink(stdout_path)
    return stdout_content

def format_linkedread(name, bc, seq, qual, forward: bool):
    '''Given a linked-read output type, will format the read accordingly and return it'''
    if c.outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * c.barcodebp}{qual}']
        if bc not in c.used_bc:
            c.used_bc[bc] = True
            c.barcodelist.write(f"{bc}\n")
    elif c.outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}:{bc} {fr}', seq, '+', qual]
        if bc not in c.used_bc:
            c.barcodelist.write(f"{bc}\n")
            c.used_bc[bc] = True
    elif c.outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        if bc not in c.used_bc:
            acbd = "".join(next(c.bc_generator))
            c.used_bc[bc] = acbd
            c.barcodelist.write(f"{bc}\n")
        else:
            acbd = c.used_bc[bc]
        read = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{acbd}', seq, '+', qual]
    elif c.outformat == "standard":
        fr = "/1" if forward else "/2"
        if bc not in c.used_bc:
            c.used_bc[bc] = bc
            c.barcodelist.write(f"{bc}\n")
        read = [f'@{name}{fr}\tBV:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif c.outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        if bc not in c.used_bc:
            stlfr_bc = "_".join([str(i) for i in next(c.bc_generator)])
            c.used_bc[bc] = stlfr_bc
            c.barcodelist.write(f"{bc}\n")
        else:
            stlfr_bc = c.used_bc[bc]
        read = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    return read

def randomlong(Par,seq_,EXPM):
    '''Length of molecules is randomly distributed according to an exponential distribution'''
    index = 0	
    lensingle = len(seq_)
    molecules = []
    for i in range(EXPM):
        start = int(RNG.uniform(low = 0,high = lensingle))
        length = int(RNG.exponential(scale = c.mollen))
        if length != 0:
            end = start + length - 1
            if end>lensingle:
                Molseq = seq_[start:lensingle]
                lengthnew = lensingle - start
                NewMol = Molecule(lengthnew,start,lensingle,index)
                molecules.append(NewMol)
            else:
                Molseq = seq_[start:end]
                NewMol = Molecule(length - 1,start,end,index)
                molecules.append(NewMol)
            index += 1
    return molecules

def deternumdroplet(molecules,molnum):
    '''Determine the number of droplets'''
    large_droplet = c.totalbarcodes
    assign_drop = []
    if molnum < 0:
        frag_drop = [abs(molnum)]
    else:
        try:
            frag_drop = RNG.poisson(molnum,large_droplet)
        except np._core._exceptions._ArrayMemoryError as e:
            mimick_errorterminate(f"The barcode and type combination yields far too many barcodes to sample. Error as reported by numpy:\n{e}")
    totalfrag = 0
    
    for i in range(large_droplet):
        totalfrag += frag_drop[i]
        if totalfrag <= len(molecules):
            assign_drop.append(frag_drop[i])
        else:
            last = len(molecules) - (totalfrag-frag_drop[i])
            assign_drop.append(last)
            break
    return assign_drop

def selectbarcode(drop,molecules,c):
    '''Select barcode to use for droplet/partition'''
    permutnum = RNG.permutation(len(molecules))
    N_droplet = len(drop)
    assigned_barcodes = set()
    droplet_container = []
    start = 0
    
    for i in range(N_droplet):
        num_molecule_per_partition = drop[i]
        index_molecule = permutnum[start:start + num_molecule_per_partition]
        totalseqlen = 0
        temp = []
        start = start + num_molecule_per_partition
        try:
            bc = next(c.barcodes) if c.barcodetype in ["10x", "tellseq"] else "".join(next(c.barcodes))
        except StopIteration:
            mimick_errorterminate('[Error] No more barcodes left for simulation. The requested parameters require more barcodes.')
        c.remainingbarcodes -= 1
        for j in range(num_molecule_per_partition):
            index = index_molecule[j]
            temp.append(index)
            molecules[index].index_droplet = i
            molecules[index].barcode = bc
            totalseqlen += molecules[index].length
        assigned_barcodes.add(bc)
        droplet_container.append(temp)
    return droplet_container, assigned_barcodes

def MolSim(processor,molecule,_fasta,w,c):
    '''Parallelize linked reads simulation'''
    for mol in molecule:
        moleculenumber = mol.seqidx + 1
        moleculedroplet = mol.index_droplet + 1
        barcodestring = mol.barcode
        chromstart = w.start + mol.start
        chromend = w.start + mol.end

        header = f'MOL:{moleculenumber}_GEM:{moleculedroplet}_BAR:{barcodestring}_CHROM:{w.chrom}_START:{chromstart}_END:{chromend}'
        try:
            seq__ = _fasta[w.chrom][w.start+mol.start-1:w.start+mol.end-1].seq
        except Exception as e:
            mimick_errorterminate(f"[Error] There was an issue with pyfaidx accessing the genomic interval requested. Error reported by pyfaidx:\n{e}", False)
        truedim = mol.length-seq__.count('N')
        if c.molcov < 1:
            N = int(truedim*c.molcov)/(c.length*2)
        else:
            # draw N from a normal distribution with a mean of molcov and stdev of molcov/3, avoiding < 0
            N = max(0, int(RNG.normal(c.molcov, c.molcov/3)))
            # set ceiling to avoid N being greater than can be sampled
            N = min(N, int(truedim/(c.length*2)))
        R1A = os.path.abspath(f'{c.OUT}/{c.PREFIX}_S1_L{c.hapnumber.zfill(3)}_R1_001.fastq')
        R2A = os.path.abspath(f'{c.OUT}/{c.PREFIX}_S1_L{c.hapnumber.zfill(3)}_R2_001.fastq')

        if N != 0:
            molfa = os.path.abspath(f'{c.OUT}/{processor}_{moleculenumber}.fa')
            with open(molfa, 'w') as faout:
                faout.write(f'>{header}\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')
            R1tmp = os.path.abspath(f"{c.OUT}/{processor}.R1.tmp.fq")
            R2 = os.path.abspath(f"{c.OUT}/{processor}.R2.fq")
            try:
                stdout = wrap_wgsim(
                    r1 = R1tmp,
                    r2 =R2,
                    ref = molfa,
                    err_rate = c.error,
                    mut_rate = c.mutation,
                    indel_frac = c.indels,
                    indel_ext = c.extindels,
                    N = N,
                    dist = c.distance,
                    stdev = c.stdev,
                    size_l = c.len_r1,
                    size_r = c.len_r2,
                    max_n = 0.05,
                    is_hap = 0,
                    is_fixed = 0,
                    seed = 0
                )
                with open(f'{c.OUT}/{c.PREFIX}.wgsim.mutations', 'a') as f:
                    f.write(stdout)
            except KeyboardInterrupt:
                mimick_keyboardterminate()

            os.remove(molfa)
            if os.stat(R1tmp).st_size == 0:
                os.remove(R1tmp)
                os.remove(R2)
            else:
                with open(R1tmp,'r') as infile, open(R1A,'a') as outfile:
                    for name,seq,qual in readfq(infile):
                        read = format_linkedread(name, barcodestring, seq, qual, True)
                        outfile.write('\n'.join(read) + '\n')
                os.remove(R1tmp)
                with open(R2,'r') as infile, open(R2A,'a') as outfile:
                    for name,seq,qual in readfq(infile):
                        if c.outformat == "10x":
                            read = [f'@{name}',seq,'+',qual]
                        else:
                            read = format_linkedread(name, barcodestring, seq, qual, False)
                        outfile.write('\n'.join(read) + '\n')
                os.remove(R2)

def LinkedSim(w,c):
    '''Perform linked-reads simulation'''
    try:
        _fasta = pyfaidx.Fasta(c.ffile)
    except ValueError as e:
        mimick_errorterminate(f'[Error] {e} in {os.path.basename(c.ffile)}')
    if w.chrom not in _fasta:
        mimick_console.log(f'[Warning] Chromosome {w.chrom} not found in {c.ffile}. Skipping.', style = "yellow")
        return
    chr_ = _fasta[w.chrom]
    seq_ = chr_[w.start-1:w.end].seq
    region = f"{w.chrom}_{w.start}_{w.end}"
    Ns = seq_.count('N') #normalize coverage on Ns

    MRPM = (c.molcov*c.mollen)/(c.length*2) if c.molcov < 1 else c.molcov
    TOTALR = round(((c.regioncoverage*(len(seq_)-Ns))/c.length)/2)
    EXPM = round(TOTALR/MRPM)

    info_table = log_table()
    info_table.add_row('Average number of paired reads per molecule', f'{round(MRPM)}')
    info_table.add_row('Reads required to get the expected coverage', f'{TOTALR}')
    info_table.add_row('Expected number of molecules', f'{EXPM}')
    # MolSet
    molecules = randomlong(c,seq_,EXPM)
    info_table.add_row('Molecules generated', f'{len(molecules)}')
    drop = deternumdroplet(molecules,c.molnum)
    info_table.add_row('Partitions molecules assigned to', f'{len(drop)}')
    info_table.add_row('Available barcodes', f'{c.remainingbarcodes}')
    droplet_container,assigned_barcodes = selectbarcode(drop,molecules,c)
    info_table.add_row('Barcodes remaining after molecule assignment', f'{c.remainingbarcodes}')
    mimick_console.log(info_table)
    chunk_size = len(molecules)/c.threads
    slices = Chunks(molecules,math.ceil(chunk_size))

    processes=[]
    error_occurred = multiprocessing.Event()
    for i,molecule in enumerate(slices,1):
        processor = f'p{i}'
        p = multiprocessing.Process(target=MolSim, args=(processor,molecule,_fasta,w,c), daemon=True)
        p.start()
        processes.append(p)

    try:
        while processes:
            for p in processes[:]:  # Iterate over a copy of the list
                if not p.is_alive():  # Process finished
                    p.join()  # Get exit code
                    if p.exitcode != 0:  # Process failed
                        error_occurred.set()
                    processes.remove(p)
            
            if error_occurred.is_set():
                mimick_console.rule("Terminating Mimick due to an error", style = "red")
                break
            
            time.sleep(0.1)  # Small delay to prevent busy waiting

    finally:
        # Cleanup - terminate any remaining processes if error occurred
        if error_occurred.is_set():
            for p in processes:
                if p.is_alive():
                    p.terminate()
        
        # Wait for all processes to finish (either normally or terminated)
        for p in processes:
            p.join()

    if error_occurred.is_set():
        sys.exit(1)