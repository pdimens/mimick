#! /usr/bin/env python3

import queue
import os
import shutil
import pysam
import threading
from .common import mimick_console, mimick_keyboardterminate
from .file_ops import readfq

stop_event = threading.Event()

def format_linkedread(name, bc, outbc, outformat, seq, qual, forward: bool):
    '''Given a linked-read output type, will format the read accordingly and return it'''
    if outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * len(bc)}{qual}']
    elif outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name}:{bc} {fr}', seq, '+', qual]
    elif outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        sequence = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{outbc}', seq, '+', qual]
    elif outformat == "standard":
        fr = "/1" if forward else "/2"
        sequence = [f'@{name}{fr}\tVX:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        sequence = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    return "\n".join(sequence)

# Worker function
def append_worker(R1_fq, R2_fq, gff, output_format, quiet, queue):
    '''
    This worker function runs on its own thread and monitors a queue
    of input files to append to the output files, then deleting the input files. It's intended to be
    initialized with the final R1/R2 and gff files (that will be appended to),
    along with the data output format ('stlfr`, tellseq`, etc.) and the queue to
    monitor. The queue is expected to contain:
    (temporary_R1.fq, temporary_R2.fq, temporary_GFF, native barcode, converted barcode)

    If the queue received a None value, it will initiate the final process of bgzipping
    R1_fq and R2_fq and deleting the uncompressed files.
    '''
    while not stop_event.is_set():
        item = queue.get()
        if item is None:
            # Exit signal received
            if quiet < 2:
                mimick_console.log(f"Compressing [blue]{os.path.basename(R1_fq)}[/]")
            pysam.tabix_compress(R1_fq, f'{R1_fq}.gz', force=True)
            os.remove(R1_fq)
            if quiet < 2:
                mimick_console.log(f"Compressing [blue]{os.path.basename(R2_fq)}[/]")
            pysam.tabix_compress(R2_fq, f'{R2_fq}.gz', force=True)
            os.remove(R2_fq)
            queue.task_done()
            break

        temp1, temp2, temp_gff, barcode, output_barcode = item
        try:
            with open(R1_fq, 'a') as out, open(temp1, 'r') as src:
                for name,seq,qual in readfq(src):
                    sequence = format_linkedread(
                        name = name,
                        bc = barcode,
                        outbc = output_barcode,
                        outformat = output_format,
                        seq = seq,
                        qual = qual,
                        forward = True
                    )
                    out.write(sequence + '\n')
            
            with open(R2_fq, 'a') as out, open(temp2, 'r') as src:
                for name,seq,qual in readfq(src):
                    if output_format == "10x":
                        sequence = '\n'.join([f'@{name}',seq,'+',qual])
                    else:
                        sequence = format_linkedread(
                            name = name,
                            bc = barcode,
                            outbc = output_barcode,
                            outformat = output_format,
                            seq = seq,
                            qual = qual,
                            forward = False
                        )
                    out.write(sequence + '\n')

            with open(temp_gff, 'r') as src, open(gff, 'a') as dest:
                shutil.copyfileobj(src, dest)

        except Exception as e:
            print(f"Error processing {temp1}, {temp2}, {temp_gff}: {e}")
            stop_event.set()
        finally:
            os.remove(temp1)
            os.remove(temp2)
            os.remove(temp_gff)
            queue.task_done()
