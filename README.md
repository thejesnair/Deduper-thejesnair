# Deduper

## Reference Based PCR Duplicate Removal Tool

`nair-deduper.py` is a Python script that removes PCR duplicates from sorted SAM files using UMIs and retains a single copy for each read. 


### Deduplication Logic
Duplicates are identified based on:
- chromosome
- strand
- adjusted genomic position
- CIGAR string
- UMI

Genomic positions are adjusted for soft clipping based on CIGAR string. Hard clipping is not addressed in this script.


### Input
- sorted, uniquely mapped SAM file
- UMI file containing valid UMIs (8 bp long)

### Output
- Deduplicated SAM file
- Summary statistics printed to terminal with chromosomes sorted as numeric, X/Y/MT, then scaffolds

### Usage
- Script requires argparse inputs:
    - `-f`, SAM input file path
    - `-o`, output path for deduplicated SAM file
    - `-u`, file containing list of UMIs
- to run script, please use the following command in the terminal
`./nair_deduper.py -u <umi.txt> -f <in.sam> -o <out.sam>`