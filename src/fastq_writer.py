#! /usr/bin/env python3

import os
import threading
import queue
from .file_ops import readfq

def format_linkedread(name, bc, outbc, outformat, seq, qual, forward: bool):
    '''Given a linked-read output type, will format the read accordingly and return it'''
    if outformat == "10x":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name} {fr}', f"{bc}{seq}", '+', f'{qual[0] * len(bc)}{qual}']
    elif outformat == "tellseq":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}:{bc} {fr}', seq, '+', qual]
    elif outformat == "haplotagging":
        fr = "/1" if forward else "/2"
        read = [f'@{name}{fr}\tOX:Z:{bc}\tBX:Z:{outbc}', seq, '+', qual]
    elif outformat == "standard":
        fr = "/1" if forward else "/2"
        read = [f'@{name}{fr}\tVX:i:1\tBX:Z:{bc}', seq, '+', qual]
    elif outformat == "stlfr":
        fr = "1:N:0:ATAGCT" if forward else "2:N:0:ATAGCT"
        read = [f'@{name}#{stlfr_bc} {fr}', seq, '+', qual]
    return read


# Create a thread-safe queue to hold (temp1, temp2) tuples
WRITER_QUEUE = queue.Queue()

# Worker function
def append_worker(R1_fq, R2_fq, output_format):
    while True:
        item = WRITER_QUEUE.get()
        if item is None:
            break  # Exit signal received

        temp1, temp2, barcode, output_barcode = item
        try:
            with open(R1_fq, 'a') as out1, open(temp1, 'r') as src:
                for name,seq,qual in readfq(src):
                    read = format_linkedread(
                        name = name,
                        bc = barcode,
                        outbc = output_barcode,
                        outformat = output_format,
                        seq = seq,
                        qual = qual,
                        forward = True
                    )
                    out1.write('\n'.join(read) + '\n')
            os.remove(temp1)
            
            with open(R2_fq, 'a') as out2, open(temp2, 'r') as src:
                for name,seq,qual in readfq(src):
                    if output_format == "10x":
                        read = [f'@{name}',seq,'+',qual]
                    else:
                        read = format_linkedread(
                            name = name,
                            bc = barcode,
                            outbc = output_barcode,
                            outformat = output_format,
                            seq = seq,
                            qual = qual,
                            forward = False
                        )
                    out2.write('\n'.join(read) + '\n')
            os.remove(temp2)   

        except Exception as e:
            print(f"Error processing {temp1}, {temp2}: {e}")
        finally:
            WRITER_QUEUE.task_done()
