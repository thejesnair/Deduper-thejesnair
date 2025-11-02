#!/usr/bin/env python

# pathway to actual deduper file: /projects/bgmp/shared/deduper
# C1_SE_uniqAlign.sam

import argparse
import gzip
import re
from collections import Counter

'''
nair_deduper.py removes PCR duplicates from a SAM file and writes sequences to an outputted SAM file

note: script is written for uniquely mapped sam file and must be sorted prior to execution
'''

def get_args():
    '''Takes in command line arguments for path to sorted input SAM file, path to output SAM, 
    path to UMI list needed for removal of PCR duplicates.'''
    parser = argparse.ArgumentParser(description=
                                     "Program to deduplicate PCR duplicates from UNIQUELY MAPPED and SORTED SAM file")
    parser.add_argument("-f", "--file", help="uniquely sorted SAM file, absolute path", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="deduped SAM file, absolute path", type=str, required=True)
    parser.add_argument("-u", "--umi", help="umi list, absolute path", type=str, required=True)

    return parser.parse_args()

args = get_args()
in_file = args.file
out_file = args.outfile
umi_file = args.umi

# cigar variables
CIGAR_REGEX = re.compile(r'(\d+)([MDNS=X])') # creates regex obj
REF_OPS = {'M', 'D', 'N', '=', 'X'} # operations that are ref consuming

# HELPER FUNCTIONS


# check known umi
def known_umi(umi, umi_set) -> bool:
    """ check if umi is in set of known umis, return T/F """
    return umi in umi_set

# extract key, create record (qname, umi, cigar)
def parse_SAM_line(line:str) -> list:
    """ Takes in SAM line string
    Returns tuple of QNAME, UMI, strand (flag 16), rname, pos, cigar string or None for header """
#for line in sam_in:
    if line.startswith('@'):
        return None # skip header
    fields = line.strip().split("\t")
    qname = fields[0]
    umi = qname.split(":")[-1]
    flag = int(fields[1])
    rname = fields[2]
    pos = int(fields[3]) # 1-based leftmost mapping position
    cigar_string = fields[5]
    return qname, umi, flag, rname, pos, cigar_string

# get int, op pairs from cigar string
def split_cigar(cigar) -> list[tuple[int, str]]:
    """ return n, op for cigar string 
    ex: "3S20M" -> [(3, 'S'), (20, 'M')] """
    return [(int(n), op) for n, op in CIGAR_REGEX.findall(cigar)]

# determine len of read
def ref_consumed_len(ops) -> int:
    """ sum only ref consuming operations, NOT soft clipping """
    return sum(n for n, op in ops if op in REF_OPS) # total len

# determine soft-clipped lens
def soft_clipping(cigar_ops) -> tuple[int, int]:
    """ determine soft clipped regions if applicable
    ex: list of (n, op) from split_cigar
    [(5, 'S'), (20, 'M')] -> (5, 0)
    [(30, 'M'), (4, 'S')] -> (0, 4) """
    S_left = 0
    S_right = 0
    for (n, op) in cigar_ops:
        if cigar_ops[0][1] == 'S': #where left hand soft clipping
            S_left = cigar_ops[0][0] #n
        if cigar_ops[-1][1] == 'S':
            S_right = cigar_ops[-1][0] #n
    return S_left, S_right

# determine strand
def assign_strand(flag) -> bool:
    """ return +/- for forward or reverse strand """
    return (flag & 16) == 16

# calculate adjusted 5' position
def adj_5prime(pos, cigar, flag):
    """ takes in flag, cigar string, and flag
    returns adjusted 5' position
    """
    S_left, S_right = soft_clipping(split_cigar(cigar))
    if assign_strand(flag) is False:
        return pos - S_left
    else:
        ref_len = ref_consumed_len(split_cigar(cigar))
        return pos + ref_len - 1 + S_right
    
# create a key to match duplicates
def record_key(rname, flag, pos, cigar, umi):
    """ creates a key of each sam record
    if keys are exact matches--> duplicate
    if one parameter differs then both records are kept """
    strand_dir = assign_strand(flag)
    adj_5pos = adj_5prime(pos, cigar, flag)
    return (rname, strand_dir, adj_5pos, umi) # error here earlier, don't include cigar bc reads can have diff cigars

# MAIN FXN
def dedup(in_file, umi_file, out_file):
    """ uses helper fxns to parse through sam, create known UMI set
    record key of reads, ignore duplicates, write sam out file """
    # create UMI set
    with open(umi_file, 'r') as umi_file:
        umi_set = {line.strip() for line in umi_file}

    current_rname = None # track when read (contig) changes, can clear seen_records
    seen_records = set()

    # COUNTERS for output report

    header_count = 0
    wrong_umi_count = 0
    unique_reads_count = 0
    removed_dupes_count = 0
    reads_per_chrom = Counter()

    # write headers to out sam file, grab first read
    with open(in_file, 'r') as sam_in, open(out_file, 'w') as sam_out:
        first_record = None
        for line in sam_in:
            if line.startswith("@"):
                header_count += 1
                sam_out.write(line)
            else:
                first_record = line
                break

    # process reads and PCR dupes
        def compare_reads(line):
            nonlocal current_rname, wrong_umi_count, unique_reads_count, removed_dupes_count
            # nonlocal to access w/i nested fxn
            record = parse_SAM_line(line)
            qname, umi, flag, rname, pos, cigar_string = record

            # if moved onto new read, clear set!
            if rname != current_rname:
                current_rname = rname
                seen_records.clear() # resets memory
            if not known_umi(umi, umi_set):
                wrong_umi_count += 1
                return
            key = record_key(rname, flag, pos, cigar_string, umi)
            if key in seen_records:
                removed_dupes_count += 1
                return
            
            seen_records.add(key)
            sam_out.write(line)
            unique_reads_count += 1
            reads_per_chrom[rname] += 1
            
        if first_record:
            compare_reads(first_record)

        for line in sam_in:
            compare_reads(line) 
    # bash command diff to see diff b/w two test output files (exp vs actual output)

    return header_count, wrong_umi_count, unique_reads_count, removed_dupes_count, reads_per_chrom 


def chrom_sort_key(chrom):
    """ sorting reads per chrom so chr 1-23, X, Y, MT, and scaffolds"""
    if chrom.isdigit():
        return int(chrom)
    if chrom == "X":
        return 20   # Mus musculus have 19 autuosomal chr, '20' would be X
    if chrom == "Y":
        return 21
    if chrom == "MT":
        return 22
    return 23  # scaffolds go last

if __name__ == "__main__":
    # print summary for qualtrics submission
    header_count, wrong_umi_count, unique_reads_count, removed_dupes_count, reads_per_chrom = dedup(in_file, umi_file, out_file)

    print(f"Number of header lines: {header_count}")
    print(f"Number of unique reads: {unique_reads_count}")
    print(f"Number of wrong UMIs: {wrong_umi_count}")
    print(f"Number of removed duplicates: {removed_dupes_count}")

    for chrom in sorted(reads_per_chrom, key=chrom_sort_key):
        print(f"{chrom}\t{reads_per_chrom[chrom]}")
    
        



