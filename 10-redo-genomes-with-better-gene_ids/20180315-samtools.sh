#!/bin/sh
## Set the job name
#PBS -N 20180315-samtools
#PBS -l nodes=1:ppn=4,vmem=30gb
# Run my job

module load CBC
module load samtools

echo "transopsonome"
echo "25 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/25transposonome
samtools view -h -f 0x2 -b 25_tspn_aligned.sortedByCoord.out.bam > 25_tspn_aligned.sortedbycoord.f2.bam
samtools view 25_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "26 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/26transposonome
samtools view -h -f 0x2 -b 26_tspn_aligned.sortedByCoord.out.bam > 26_tspn_aligned.sortedbycoord.f2.bam
samtools view 26_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "27 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/27transposonome
samtools view -h -f 0x2 -b 27_tspn_aligned.sortedByCoord.out.bam > 27_tspn_aligned.sortedbycoord.f2.bam
samtools view 27_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "28 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/28transposonome
samtools view -h -f 0x2 -b 28_tspn_aligned.sortedByCoord.out.bam > 28_tspn_aligned.sortedbycoord.f2.bam
samtools view 28_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "29 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/29transposonome
samtools view -h -f 0x2 -b 29_tspn_aligned.sortedByCoord.out.bam > 29_tspn_aligned.sortedbycoord.f2.bam
samtools view 29_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "30 transposonome"
cd /home/sf040090/Li-2017-hFGC/cln_tspn_19W/30transposonome
samtools view -h -f 0x2 -b 30_tspn_aligned.sortedByCoord.out.bam > 30_tspn_aligned.sortedbycoord.f2.bam
samtools view 30_tspn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "transcriptome"
echo "25 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/25transcriptome
samtools view 25_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 25_txn_aligned.sortedbycoord.out.bam > 25_txn_aligned.sortedbycoord.f2.bam
samtools view 25_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "26 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/26transcriptome
mv Aligned.sortedByCoord.out.bam 26_txn_aligned.sortedbycoord.out.bam
samtools view 26_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 26_txn_aligned.sortedbycoord.out.bam > 26_txn_aligned.sortedbycoord.f2.bam
samtools view 26_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "27 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/27transcriptome
mv Aligned.sortedByCoord.out.bam 27_txn_aligned.sortedbycoord.out.bam
samtools view 27_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 27_txn_aligned.sortedbycoord.out.bam > 27_txn_aligned.sortedbycoord.f2.bam
samtools view 27_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "28 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/28transcriptome
mv Aligned.sortedByCoord.out.bam 28_txn_aligned.sortedbycoord.out.bam
samtools view 28_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 28_txn_aligned.sortedbycoord.out.bam > 28_txn_aligned.sortedbycoord.f2.bam
samtools view 28_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "29 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/29transcriptome
mv Aligned.sortedByCoord.out.bam 29_txn_aligned.sortedbycoord.out.bam
samtools view 29_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 29_txn_aligned.sortedbycoord.out.bam > 29_txn_aligned.sortedbycoord.f2.bam
samtools view 29_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "30 transcriptome"
cd /home/sf040090/Li-2017-hFGC/cln_txn_19W/30transcriptome
mv Aligned.sortedByCoord.out.bam 30_txn_aligned.sortedbycoord.out.bam
samtools view 30_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c
samtools view -h -f 0x2 -b 30_txn_aligned.sortedbycoord.out.bam > 30_txn_aligned.sortedbycoord.f2.bam
samtools view 30_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c

echo "done"