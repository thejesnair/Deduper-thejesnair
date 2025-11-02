#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=run_dedup
#SBATCH --output=run_dedup%j

python_script="/projects/bgmp/tnair/bioinfo/Bi624/Deduper-thejesnair/nair_deduper.py"
sam_file="/projects/bgmp/tnair/bioinfo/Bi624/Deduper-thejesnair/C1_SE_uniqAlign.sam"
in_sorted_sam="/projects/bgmp/tnair/bioinfo/Bi624/Deduper-thejesnair/C1_SE_uniqAlign_sorted.sam"
umi_file="/projects/bgmp/tnair/bioinfo/Bi624/Deduper-thejesnair/STL96.txt"
out_sam="/projects/bgmp/tnair/bioinfo/Bi624/Deduper-thejesnair/C1_SE_uniqAlign_sorted_deduped.sam"

# activate conda envs, sort, deactivate 
echo "Activating samtools env"
conda activate samtools

echo "Sorting SAM file"
/usr/bin/time -v samtools sort $sam_file -o $in_sorted_sam

conda deactivate

# dedupe
/usr/bin/time -v python $python_script -f $in_sorted_sam -u $umi_file -o $out_sam