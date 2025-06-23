# Mimick
Linked-read sequence simulator

![mimick_logo](_media/mimick_logo.png)

Mimick is a simulator for linked-read FASTQ data.Mimick started its life as XENIA, a module in
[VISOR](https://github.com/davidebolo1993/VISOR).  It allows you to simulate an arbitrary number
of haplotypes, set overall coverage, molecule coverage, and mix-match barcodes with linked-read chemistries.

## Supported Linked-Read Types:
- 10X
- Haplotagging
- stLFR
- TELLseq

## Get Started
See the [installation guide](install.md) and then call up `mimick` in the command line to be greeted with:

```terminal
$| mimick

 Usage: mimick [OPTIONS] BARCODES FASTA...                                                              

 Simulate linked-read FASTQ using genome haplotypes. Barcodes can be supplied one of two ways:          

  1 let Mimick randomly generate barcodes based on a specification of length,count                      
     • two integers, comma-separated, no space                                                          
     • e.g. 16,400000 would generate 400,000 unique 16bp barcodes                                       
  2 you can provide a file of specific nucleotide barcodes, 1 per line                                  

 You can specify the linked-read barcode chemistry to simulate via --lr-type as well as the output      
 format of FASTQ files (default is the same as barcode type). For example, you can generate 96 barcodes 
 (common haplotagging style), select --lr-type stlfr (combinatorial 3-barcode on R2 read), and have     
 --output-type tellseq (@seqid:barcode header format).                                                  


   --lr-type      Format                                                                                
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━                               
   10x/tellseq    single barcode on R1                                                                  
   haplotagging   R1 and R2 each have different combinatorial 2-barcodes                                
   stlfr          combinatorial 3-barcode on R2                                                         


   --output-type   Barcode Location                                      Example                        
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━      
   10x             start of R1 sequence                                  ATAGACCATAGAGGACA...           
   haplotagging    sequence header as BX:Z:ACBD                          @SEQID BX:Z:A0C331B34D87       
   standard        sequence header as BX:Z:BARCODE, no specific format   @SEQID BX:Z:ATACGAGACA         
   stlfr           appended to sequence ID via #1_2_3                    @SEQID#1_354_39                
   tellseq         appended to sequence ID via :ATCG                     @SEQID:TATTAGCAC               


╭─ General Options ────────────────────────────────────────────────────────────────────────────────────╮
│ --help               Show this message and exit.                                                     │
│ --circular       -C  contigs are circular/prokaryotic                                                │
│ --output-prefix  -o  output file prefix                                                              │
│                      [default: simulated/SIM]                                                        │
│ --output-type    -O  output format of FASTQ files                                                    │
│ --quiet          -q  0 all output, 1 no progress bar, 2 no output                                    │
│                      [default: 0]                                                                    │
│ --regions        -r  one or more regions to simulate, in BED format                                  │
│ --seed           -S  random seed for simulation                                                      │
│ --threads        -t  number of threads to use for simulation                                         │
│                      [default: 2; x>=1]                                                              │
│ --version            Show the version and exit.                                                      │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Read Simulation Parameters ─────────────────────────────────────────────────────────────────────────╮
│ --coverage     mean coverage target for simulated data                                               │
│                [default: 30.0; x>=0.05]                                                              │
│ --distance     outer distance between the two ends in bp                                             │
│                [default: 500; x>=0]                                                                  │
│ --error        base error rate                                                                       │
│                [default: 0.02; 0<=x<=1]                                                              │
│ --extindels    indels extension rate                                                                 │
│                [default: 0.25; 0<=x<=1]                                                              │
│ --indels       indels rate                                                                           │
│                [default: 0.15; 0<=x<=1]                                                              │
│ --lengths      length of R1,R2 reads in bp                                                           │
│                [default: 150,150]                                                                    │
│ --mutation     mutation rate                                                                         │
│                [default: 0.001; x>=0]                                                                │
│ --stdev        standard deviation of --distance                                                      │
│                [default: 50; x>=0]                                                                   │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Linked Read Parameters ─────────────────────────────────────────────────────────────────────────────╮
│ --combinatorial      -b  treat barcodes as combinatorial with this many segments                     │
│                          [default: haplotagging]                                                     │
│ --molecule-attempts  -a  how many tries to create a molecule with <70% ambiguous bases               │
│                          [default: 300; x>=5]                                                        │
│ --molecule-coverage  -c  mean percent coverage per molecule if <1, else mean number of reads per     │
│                          molecule                                                                    │
│                          [default: 0.2; x>=1e-05]                                                    │
│ --molecule-length    -m  mean length of molecules in bp                                              │
│                          [default: 80000; x>=650]                                                    │
│ --molecules-per      -n  mean number of unrelated molecules per barcode per chromosome, where a      │
│                          negative number (e.g. -2) will use a fixed number of unrelated molecules    │
│                          and a positive one will draw from a Normal distribution                     │
│                          [default: 2]                                                                │
│ --singletons         -s  proportion of barcodes will only have a single read pair                    │
│                          [default: 0; 0<=x<=1]                                                       │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────╯


 Documentation: https://pdimens.github.io/mimick/
```

!> Why "mimick"? It's simple, really. This software mimics linked-read data, Pavel has an affinity for naming software after [fictional monsters](https://en.wikipedia.org/wiki/Mimic_(Dungeons_%26_Dragons)) and "mimick" (with a "k") is the old-English spelling of the word, leaving `mimic` available for some other bioinformatician to use for a less farcical reason. Despite the lore of mimics being deadly traps, this software is anything but, we promise.

### Authors
<img src="https://avatars.githubusercontent.com/u/19176506?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@pdimens](https://github.com/pdimens) (Mimick)

<img src="https://avatars.githubusercontent.com/u/39052119?v=4" width="50" height="50" style="border-radius: 50%; object-fit: cover;"/> [@davidebolo1993](https://github.com/davidebolo1993) (VISOR)



