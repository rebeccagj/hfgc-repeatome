#!/bin/sh
## Set the job name
#PBS -N sort_index_geneinfo_19W
#PBS -l nodes=1:ppn=1,vmem=30gb
# Run my job

mkdir /scratch/sf040090/indexsort
cd /scratch/sf040090/indexsort

samtools sort /home/sf040090/Li-2017-hFGC/mapped/25_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 25_assigned_sorted.bam
samtools index 25_assigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/mapped/26_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 26_assigned_sorted.bam
samtools index 26_assigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/mapped/27_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 27_assigned_sorted.bam
samtools index 27_assigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/mapped/28_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 28_assigned_sorted.bam
samtools index 28_assigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/mapped/29_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 29_assigned_sorted.bam
samtools index 29_assigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/mapped/30_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 30_assigned_sorted.bam
samtools index 30_assigned_sorted.bam

mv /scratch/sf040090/indexsort /home/sf040090/Li-2017-hFGC/mapped

echo "onward to step 6!"