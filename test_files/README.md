This folder contains test files used to develop `nair-deduper.py`

## Files
- test_sam_input.sam
    - contains test reads, formatted like a SAM file but simplified for readability and quick testing (SORTED)
- test_sam_output.sam
    - expected output file generated from input
    - duplicates should be removed and the original records kept
- test_umi.txt
    - file containing known UMIs


## Cases being tested
1. Unique, no duplicates
    - `readU_xxx:CCGTTAAC`
    - `read6_unique_xxx:CCGTTAAC`
    - `readMT_unique_xxx:ATCCATGG`
2. Forward read soft clip: SAME adj 5' pos, SAME UMI (duplicate, keep first instance)
    - `readA1_fwd_xxx:GAACAGGT`
    - `readA2_fwd_xxx:GAACAGGT`
    - `read6_dup1_fwd_xxx:GAACAGGT`
    - `read6_dup2_fwd_xxx:GAACAGGT`
3. Reverse read soft clip: SAME adj 5' pos, SAME UMI (duplicate, keep first instance)
    - `readC1_rev_xxx:GAACAGGT`
    - `readC2_rev_xxx:GAACAGGT`
4. SAME chr (site), SAME strand (+ or -), DIFF UMI (not duplicate, keep both)
    - `readD2_xxx:TTTACGGA`
    - `readD1_xxx:GAACAGGT`
5. SAME UMI, SAME chr, DIFF strand (not duplicate, keep both)
    - `readK_rev_xxx:ATCCATGG`
    - `readK_fwd_xxx:ATCCATGG`
6. Invalid UMI (not in known list) (discard)
    - `readX_xxx:INVALIDU`
    - `readMT_badumi_xxx:INVALIDU`
