#!/bin/sh
## Set the job name
#PBS -N 20180315-txn-samtools
#PBS -l nodes=1:ppn=4,vmem=30gb
# Run my job

module load CBC
module load samtools

echo "25 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/25transcriptome
mv Aligned.sortedByCoord.out.bam 25_txn_aligned.sortedByCoord.out.bam
samtools view -h -f 0x2 25_txn_aligned.sortedByCoord.out.bam > 25_txn_aligned.sortedbycoord.f2.bam
# note that multi-mappers are kept!!!!! for non-transposon data might be good to remove multi-mappers (e.g., -Fh 100)

echo "26 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/26transcriptome
mv Aligned.sortedByCoord.out.bam 26_txn_aligned.sortedByCoord.out.bam
samtools view 26_txn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 26_txn_aligned.sortedByCoord.out.bam > 26_txn_aligned.sortedbycoord.f2.bam

echo "27 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/27transcriptome
mv Aligned.sortedByCoord.out.bam 27_txn_aligned.sortedByCoord.out.bam
samtools view 27_txn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 27_txn_aligned.sortedByCoord.out.bam > 27_txn_aligned.sortedbycoord.f2.bam

echo "28 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/28transcriptome
mv Aligned.sortedByCoord.out.bam 28_txn_aligned.sortedByCoord.out.bam
samtools view 28_txn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 28_txn_aligned.sortedByCoord.out.bam > 28_txn_aligned.sortedbycoord.f2.bam

echo "29 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/29transcriptome
mv Aligned.sortedByCoord.out.bam 29_txn_aligned.sortedByCoord.out.bam
samtools view 29_txn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 29_txn_aligned.sortedByCoord.out.bam > 29_txn_aligned.sortedbycoord.f2.bam

echo "30 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/30transcriptome
mv Aligned.sortedByCoord.out.bam 30_txn_aligned.sortedByCoord.out.bam
samtools view 30_txn_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 30_txn_aligned.sortedByCoord.out.bam > 30_txn_aligned.sortedbycoord.f2.bam