#!/bin/sh
## Set the job name
#PBS -N count-molecules-wide2
#PBS -l nodes=1:ppn=8,vmem=30gb
# Run my job

mkdir /scratch/sf040090/UMIcountwide2
cd /scratch/sf040090/UMIcountwide2

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/25_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/25_counts.tsv.gz

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/26_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/26_counts.tsv.gz

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/27_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/27_counts.tsv.gz

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/28_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/28_counts.tsv.gz

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/29_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/29_counts.tsv.gz

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I /home/sf040090/Li-2017-hFGC/mapped/indexsort/30_assigned_sorted.bam -S /scratch/sf040090/UMIcountwide2/30_counts.tsv.gz

mv /scratch/sf040090/UMIcountwide2 /home/sf040090/Li-2017-hFGC/mapped

echo "donedone?"