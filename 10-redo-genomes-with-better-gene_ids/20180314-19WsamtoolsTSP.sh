#!/bin/sh
## Set the job name
#PBS -N 20180314-19WsamtoolsTSP
#PBS -l nodes=1:ppn=8,vmem=30gb
# Run my job

module load CBC
module load samtools

# TRANSPOSONOME
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/25transposonome
mv Aligned.sortedByCoord.out.bam 25_tspn_aligned.sortedByCoord.out.bam
samtools view 25_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 25_tspn_aligned.sortedByCoord.out.bam > 25_tspn_aligned.sortedbycoord.f2.bam
# note that multi-mappers are kept!!!!! for non-transposon data might be good to remove multi-mappers (e.g., -Fh 100)

cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/26transposonome
mv Aligned.sortedByCoord.out.bam 26_tspn_aligned.sortedByCoord.out.bam
samtools view 26_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 26_tspn_aligned.sortedByCoord.out.bam > 26_tspn_aligned.sortedbycoord.f2.bam

cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/27transposonome
mv Aligned.sortedByCoord.out.bam 27_tspn_aligned.sortedByCoord.out.bam
samtools view 27_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 27_tspn_aligned.sortedByCoord.out.bam > 27_tspn_aligned.sortedbycoord.f2.bam

cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/28transposonome
mv Aligned.sortedByCoord.out.bam 28_tspn_aligned.sortedByCoord.out.bam
samtools view 28_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 28_tspn_aligned.sortedByCoord.out.bam > 28_tspn_aligned.sortedbycoord.f2.bam

cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/29transposonome
mv Aligned.sortedByCoord.out.bam 29_tspn_aligned.sortedByCoord.out.bam
samtools view 29_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 29_tspn_aligned.sortedByCoord.out.bam > 29_tspn_aligned.sortedbycoord.f2.bam

cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/30transposonome
mv Aligned.sortedByCoord.out.bam 30_tspn_aligned.sortedByCoord.out.bam
samtools view 30_tspn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -hf 2 30_tspn_aligned.sortedByCoord.out.bam > 30_tspn_aligned.sortedbycoord.f2.bam