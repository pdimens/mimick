---
label: Home
hidden: true
icon: home
---

![](static/mimick_logo.png)

Mimick is a simulator for linked-read FASTQ data. It started its life as XENIA, a module in
[VISOR](https://github.com/davidebolo1993/VISOR), and has since added features and perks. Mimick
allows you to simulate an arbitrary number of haplotypes, set overall coverage, molecule coverage,
and mix-match barcodes with linked-read chemistries. You can also supply one fasta file and a vcf
file to simulate haplotype-resolved reads for all the samples in the vcf.

==- Why "mimick"? 
It's simple, really. This software mimics linked-read data, Pavel has an affinity for naming software after [fictional monsters](https://en.wikipedia.org/wiki/Mimic_(Dungeons_%26_Dragons)) and "mimick" (with a "k") is the old-English spelling of the word, leaving `mimic` available for some other bioinformatician to use for a less farcical reason. Despite the lore of mimics being deadly traps, this software is anything but, we promise.
===


## Supported Linked-Read Types:
- 10X
- Haplotagging
- stLFR
- TELLseq

## Get Started
See the [installation](install.md) and [usage](usage.md) guides, then call up `mimick` in the command line to be greeted with:

```bash
Usage: mimick [OPTIONS] FASTA...                                                
                                                                                
Simulate linked-read FASTQ data for one or many individuals                     
There are two modes of operation:                                               
                                                                                
 1 Input one or more FASTA files (haplotypes) to simulate linked reads for a    
   single individual.                                                           
 2 Input one FASTA and VCF file to simulate linked reads for all samples in the 
   VCF file with haplotypes reflective of their SNPs and indels.                
                                                                                
With the exception of 10x, all other formats are demultiplexed. Output can be in
standard:chemistry format (e.g. standard:stlfr outputs standard format with     
stLFR-style barcodes), where the barcode is encoded as a BX:Z: tag and a VX:i   
validation tag. Below are the common linked-read chemistries (to be used in     
--format) and their configurations:                                             
                                                                                
                                                                                
  chemistry      --read-lengths   Description            FASTQ format           
 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 
  10x               134,150       single barcode on R1   barcode inline on R1   
  haplotagging      150,150       2-barcodes in I1/I2    BX:Z:AxxCxxBxxDxx tag  
  stlfr             150,108       3-barcode on R2        @SEQID#1_2_3           
  tellseq           132,150       single barcode on R1   @SEQID:ATGC            
                                                                                
                                                                                
General Options:                                                                
  --circular       -c  contigs are circular/prokaryotic                         
  --output-prefix  -o  output file prefix                                       
                       [default=simulated/]                                     
  --quiet          -q  toggle to hide progress bar                              
  --threads        -t  number of threads to use for multi-sample simulation     
                       [default=2]                                              
  --seed           -s  random seed for simulation                               
  --vcf            -v  VCF-formatted file containing genotypes from which to    
                       create per-sample haplotypes                             
  --version            Show the version and exit.                               
  --help               Show this message and exit.                              
                                                                                
Linked-Read Simulation:                                                         
  --format             -f  FASTQ output format                                  
                           [default=standard:haplotagging]                      
  --genomic-coverage   -g  mean coverage (depth) target for simulated data      
                           [default=30.0]                                       
  --insert-size        -i  outer distance between the two read ends in bp       
                           [default=500]                                        
  --insert-stdev       -d  standard deviation for --insert-size                 
                           [default=50]                                         
  --molecule-attempts  -A  how many tries to create a molecule with <70% ambiguous bases                                      
                           [default=300]                                        
  --molecule-coverage  -C  mean percent coverage per molecule if <1, else mean number of reads per molecule                         
                           [default=0.2]                                        
  --molecule-length    -L  mean length of molecules in bp                       
                           [default=80000]                                      
  --molecules-per      -N  mean number of unrelated molecules per barcode per chromosome, where a negative
                           number (e.g. -2) will use a fixed number of unrelated molecules and a positive
                           one will draw from a distribution           
                           [default=2]                                          
  --read-lengths       -l  length of R1,R2 sequences in bp                      
                           [default=150,150]                                    
  --singletons         -S  proportion of barcodes that will only have one read pair                                                 
                           [default=0]         

 Documentation: https://pdimens.github.io/mimick/
```

### Authors
<img src="https://avatars.githubusercontent.com/u/19176506?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@pdimens](https://github.com/pdimens) (Mimick)

<img src="https://avatars.githubusercontent.com/u/39052119?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@davidebolo1993](https://github.com/davidebolo1993) (VISOR)



