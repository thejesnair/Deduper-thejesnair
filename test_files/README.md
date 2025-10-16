This folder contains files to test my deduplication.py script

*Files*<br>
- test_sam_input.sam
    - contains test reads, formatted like a SAM file but simplified for readability and quick testing (SORTED)
- test_sam_output.sam
    - expected output file generated from input
    - duplicates should be removed and the original records kept
- test_umi.txt
    - file containing known UMIs


*Cases being tested:*
1. Unique, no duplicates
2. Forward read soft clip, SAME adj 5' pos, SAME UMI (duplicate, keep first instance)
3. Reverse read soft clip, SAME adj 5' pos, SAME UMI (duplicate, keep first instance)
4. SAME chr (site), SAME strand (+ or -), DIFF UMI (not duplicate, keep both)
5. SAME UMI, SAME chr, DIFF strand (not duplicate, keep both)
6. Invalid UMI (not in known list) (discard)
7. Unmapped (flag 4) (discard)
