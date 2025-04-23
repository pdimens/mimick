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
import pysam
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
    if Container.outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * Container.barcodebp}{qual}']
        if bc not in Container.used_bc:
            Container.used_bc[bc] = True
            Container.barcodelist.write(f"{bc}\n")
    elif Container.outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}:{bc} {fr}', seq, '+', qual]
        if bc not in Container.used_bc:
            Container.barcodelist.write(f"{bc}\n")
            Container.used_bc[bc] = True
    elif Container.outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        if bc not in Container.used_bc:
            acbd = "".join(next(Container.bc_generator))
            Container.used_bc[bc] = acbd
            Container.barcodelist.write(f"{bc}\n")
        else:
            acbd = Container.used_bc[bc]
        read = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{acbd}', seq, '+', qual]
    elif Container.outformat == "standard":
        fr = "/1" if forward else "/2"
        if bc not in Container.used_bc:
            Container.used_bc[bc] = bc
            Container.barcodelist.write(f"{bc}\n")
        read = [f'@{name}{fr}\tVX:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif Container.outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        if bc not in Container.used_bc:
            stlfr_bc = "_".join([str(i) for i in next(Container.bc_generator)])
            Container.used_bc[bc] = stlfr_bc
            Container.barcodelist.write(f"{bc}\n")
        else:
            stlfr_bc = Container.used_bc[bc]
        read = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    return read

def randomlong(Par,seq_,EXPM):
    '''Length of molecules is randomly distributed according to an exponential distribution'''
    index = 0	
    lensingle = len(seq_)
    molecules = []
    for i in range(EXPM):
        start = int(RNG.uniform(low = 0,high = lensingle))
        length = int(RNG.exponential(scale = Container.mollen))
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
    large_droplet = Container.totalbarcodes
    assign_drop = []
    totalfrag = 0
    
    for i in range(large_droplet):
        frag = abs(molnum) if molnum < 0 else RNG.poisson(molnum,1)[0]
        totalfrag += frag
        if totalfrag <= len(molecules):
            assign_drop.append(frag)
        else:
            last = len(molecules) - (totalfrag-frag)
            assign_drop.append(last)
            break
    return assign_drop

def selectbarcode(drop,molecules,container):
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
            bc = "".join(next(container.barcodes))
        except StopIteration:
            error_terminate('No more barcodes left for simulation. The requested parameters require more barcodes.')
        container.remainingbarcodes -= 1
        for j in range(num_molecule_per_partition):
            index = index_molecule[j]
            temp.append(index)
            molecules[index].index_droplet = i
            molecules[index].barcode = bc
            totalseqlen += molecules[index].length
        assigned_barcodes.add(bc)
        droplet_container.append(temp)
    return droplet_container, assigned_barcodes

def MolSim(processor: int, molecule: Molecule, w: Molecule, container: Container) ->None:
    '''Parallelize linked reads simulation'''
    with pysam.FastaFile(container.ffile) as _fasta:
        for mol in molecule:
            moleculenumber = mol.seqidx + 1
            moleculedroplet = mol.index_droplet + 1
            barcodestring = mol.barcode
            chromstart = w.start + mol.start
            chromend = w.start + mol.end

            header = f'MOL:{moleculenumber}_GEM:{moleculedroplet}_BAR:{barcodestring}_CHROM:{w.chrom}_START:{chromstart}_END:{chromend}'
            try:
                seq__ = _fasta.fetch(w.chrom,chromstart-1, chromend+1)
            except Exception as e:
                error_terminate(f"There was an issue with accessing {w.chrom}:{chromstart}-{chromend} requested via indexed random-access. Error reported by pysam:\n{e}\n\n[yellow]If this file was compressed, make sure it, the solution might be to decompress it, sanitize it with `seqtk seq file.fa | bgzip > file.fa.gz` and try again.", False)
            truedim = mol.length-seq__.count('N')
            if container.molcov < 1:
                N = int(truedim*container.molcov)/(container.length*2)
            else:
                # draw N from a normal distribution with a mean of molcov and stdev of molcov/3, avoiding < 0
                N = max(0, int(RNG.normal(container.molcov, container.molcov/3)))
                # set ceiling to avoid N being greater than can be sampled
                N = min(N, int(truedim/(container.length*2)))

            if N != 0:
                molfa = os.path.abspath(f'{container.OUT}/p{processor}_{moleculenumber}.fa')
                with open(molfa, 'w') as faout:
                    faout.write(f'>{header}\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')
                R1tmp = os.path.abspath(f"{container.OUT}/p{processor}.R1.tmp.fastq")
                R2tmp = os.path.abspath(f"{container.OUT}/p{processor}.R2.fastq")
                try:
                    stdout = wrap_wgsim(
                        r1 = R1tmp,
                        r2 =R2tmp,
                        ref = molfa,
                        err_rate = container.error,
                        mut_rate = container.mutation,
                        indel_frac = container.indels,
                        indel_ext = container.extindels,
                        N = N,
                        dist = container.distance,
                        stdev = container.stdev,
                        size_l = container.len_r1,
                        size_r = container.len_r2,
                        max_n = 0.05,
                        is_hap = 0,
                        is_fixed = 0,
                        seed = 0
                    )
                    with open(f'{container.OUT}/{container.PREFIX}.p{processor}.wgsim.mutations', 'a') as f:
                        f.write(stdout)
                except KeyboardInterrupt:
                    mimick_keyboardterminate()

                os.remove(molfa)
                if os.stat(R1tmp).st_size == 0:
                    os.remove(R1tmp)
                    os.remove(R2tmp)
                else:
                    R1A = os.path.abspath(f'{container.OUT}/{container.PREFIX}.hap_{container.hapnumber.zfill(3)}.R1.fq')
                    R2A = os.path.abspath(f'{container.OUT}/{container.PREFIX}.hap_{container.hapnumber.zfill(3)}.R2.fq')

                    with open(R1tmp,'r') as infile, open(R1A,'a') as outfile:
                        for name,seq,qual in readfq(infile):
                            read = format_linkedread(name, barcodestring, seq, qual, True)
                            outfile.write('\n'.join(read) + '\n')
                    os.remove(R1tmp)
                    with open(R2,'r') as infile, open(R2A,'a') as outfile:
                        for name,seq,qual in readfq(infile):
                            if container.outformat == "10x":
                                read = [f'@{name}',seq,'+',qual]
                            else:
                                read = format_linkedread(name, barcodestring, seq, qual, False)
                            outfile.write('\n'.join(read) + '\n')
                    os.remove(R2tmp)

def LinkedSim(w: Molecule, container: Container) -> None:
    '''Perform linked-reads simulation'''
    with pysam.FastaFile(container.ffile) as _fasta:
        if w.chrom not in _fasta:
            container.CONSOLE.log(f'[Warning] Chromosome {w.chrom} not found in {container.ffile}. Skipping.', style = "yellow")
            return
        seq_ = _fasta.fetch(w.chrom, w.start-1, w.end + 1)
    Ns = seq_.count('N') #normalize coverage on Ns
    MRPM = (container.molcov*container.mollen)/(container.length*2) if container.molcov < 1 else container.molcov
    MRPM = max(1,MRPM)
    TOTALR = round(((container.regioncoverage*(len(seq_)-Ns))/container.length)/2)
    EXPM = round(TOTALR/MRPM)

    info_table = log_table()
    info_table.add_row('Average number of paired reads per molecule', f'{round(MRPM)}')
    info_table.add_row('Reads required to get the expected coverage', f'{TOTALR}')
    info_table.add_row('Expected number of molecules', f'{EXPM}')
    # MolSet
    molecules = randomlong(container, seq_,EXPM)
    info_table.add_row('Molecules generated', f'{len(molecules)}')
    drop = deternumdroplet(molecules, container.molnum)
    info_table.add_row('Partitions molecules assigned to', f'{len(drop)}')
    info_table.add_row('Available barcodes', f'{container.remainingbarcodes}')
    droplet_container,assigned_barcodes = selectbarcode(drop,molecules, container)
    info_table.add_row('Barcodes remaining after molecule assignment', f'{container.remainingbarcodes}')
    container.CONSOLE.log(info_table)
    chunk_size = len(molecules)/container.threads
    slices = Chunks(molecules,math.ceil(chunk_size))

    processes=[]
    error_occurred = multiprocessing.Event()
    for i,molecule in enumerate(slices,1):
        p = multiprocessing.Process(target=MolSim, args=(i,molecule, w, container), daemon=True)
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
                container.CONSOLE.rule("Terminating Mimick due to an error", style = "red")
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