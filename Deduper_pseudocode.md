## PCR Deduplication Pseudocode

Problem:

PCR is necessary for amplifying fragments when sequencing libraries, however, the introduction of duplicates can affect our analysis. We want to remove these duplicate artifacts before further downstream analyses. 


**Initial**

1. load input sam file, umi file, create output sam<br>
1b. assign all UMIs to a set<br>
1c. write all header information to output sam<br>

**Small functions**<br>

2. extract qname from SAM, return UMI for checking
3. check if unmapped (flag & 4)
4. assign strand (flag & 16)
5. parse cigar and create list containing n, op
6. determine soft clipping (left and right)
7. adjust 5' starting position based on cigar string and soft clipping
8. create record key to keep track of duplicates, containing rname, flag, adjusted five prime position, and umi

**Main function**

9. read in input sam record, umi list, and output sam files
10. parse record and assign record values to variables: qname, flag, rname, pos, cigar,
11. use previously defined fxns to determine if mapped, if known umi
12. if passes previous step then create a record key for read.<br>
12b. if come across a duplicate record (has to be the exact same record) then drop
13. write ORIGINAL record to output sam file

**Cases**

1. No duplicate
    - not a special case, will be kept
2. same strand (both + or -), same UMI<BR>
    - different soft clipping -> DUPE<br>
    code will account for this b/c it checks strandedness, UMI, and adj 5' position
3. same strand, same 5' position, different UMI --> BOTH KEPT
    - these are different records b/c of the different UMI. code will account for this
4. different strand, same position and UMI --> BOTH KEPT
    - code accounts for this
5. invalid UMI
    - if umi not in original file, record is discarded
6. read unmapped
    - discarded



*Input files:*

- UMI text file (contains 96 UMIs)
- Input SAM (sorted by POS, increasing)

*Output file:*

- Deduped SAM

bash
```
samtools
sort input sam file if not already sorted
```

python
```python
 
"""
Note: some fxns may be combined later on
for now I have listed them separately for better understanding/tracking
"""


GLOBAL variables
seen_records: set #RNAME, adj 5' pos, strand, UMI
umi_list : set of known umis

ARGPARASE
function argparse:
""" takes in filename/pathway for input sam, output sam, UMI list """

SETUP
fxn open_files(sam_input, umi, sam_output):
    open in_sam for reading:
    open UMI_list
        read each line
        put UMIs in set for set checking later
    open out_sam for writing:
        write all header information (@HD, @SQ...)

SMALL FXNs FOR READABILITY (will incorporate in main dedup fxn)
fxn parse_UMI(qname):
"""reads qname and returns UMI found at the end
ex:
input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC
return: CTGTTCAC
"""
    return UMI

fxn is_unmapped(flag):
""" returns T/F if read is unmapped """
    return (flag & 4) == 4: #unmapped, will discard

fxn assign_strand(flag):
""" returns - or + for forward or reverse strand """
    strand = #strandedness +/-
    if (flag & 16) == 16: #rev comp = TRUE
        strand = '-'
    else:
        strand = '+'
    return strand

fxn parse_cigar(cigar):
    """ return n, op for cigar string 
        ex: "3S20M" -> [(3, 'S'), (20, 'M')] """
    return list of len and char

fxn cigar_string_len(cigar):
""" returns total len of read based on cigar string input 
ex: "10M1I5M2D10M"
return 27
"""
    total_len: #len of read to determine left based pos
    for (n, op) in parse_cigar(cigar): #int and op 
        if op is in ['M', 'D', '=', 'X', 'N']:
        #cases where op consumes reference bases ^
            total_len += n
    return total_len

fxn soft_clipped(cigar):
    """ determine soft clipped regions if applicable 
    ex: [(5, 'S'), (20, 'M')]
    S_left: return 5
    [(30, 'M'), (4, 'S')]
    S_right: return 4
    """
    S_left =
    S_right = 
    cigar_list = parse_cigar(cigar)
    for (n, op) in cigar_list: #where S is first in str
        if op = S and n > 0:
            S_left = (n)
    for (n, op) in cigar_list: #where S is last in str
        if op= S and n > 0:
            S_right = (n)
    return S_left, S_right


fxn adjusted_5prime(pos, cigar, flag):
"""fxn takes in POS (col 4), CIGAR string, and FLAG(2)"""
    S_left, S_right = soft_clipped(cigar)
    if (flag & 16) != 16:
        return pos - S_left
    else:
        return pos + cigar_string_len(cigar) - 1 + S_right #subtract 1 bc already starting from position 1

fxn record_key(rname, flag, pos, cigar, umi):
"""creates a key of each record which will be used to determine duplicates. if keys are exact matches:duplicate. if one of the parameters differs then both records will be kept later"""
    strand = assign_strand(flag)
    five_prim_pos = adjusted_5prime(pos, cigar, flag)
    return (rname, strand, five_prim_pos, umi)

MAIN DEDUP FXN
fxn read SAM record(input sam, umi list, output sam):
    fxn open_files(input sam, umi list, output sam)
    seen_records = set(_)

    for line in sam line:
        if header:
            skip, already wrote to file
        parse record and split to grab each col and assign to fields:
            qname =
            flag = 
            rname = 
            pos = 
            cigar =
            already grabbed UMI previously

    #filter and discard before identifying dupes
        if is_unmapped():
            continue
        UMI = parse_umi(qname)
        if umi not in known umi list:
            continue
    #create record, THIS IS WHERE WE CHECK FOR DUPES
        key = record_key(RNAME, flag, strand, adjusted_5prime, UMI) #will filter out dupes b/c utilizing ADJUSTED 5' position and umi as well
        if key in seen_records:
            continue #drop dupe
        else:
            add key to seen_records
            write record to out_sam file #write entire original record to sam file
            (will look like out_sam.write(line))
```



