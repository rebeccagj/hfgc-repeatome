 #!/bin/sh
## Set the job name
#PBS -N 20180317-sort_index_featurecounted_BAM
#PBS -l nodes=1:ppn=4,vmem=30gb
# Run my job

module load CBC
module load samtools

mkdir /scratch/sf040090/genebycell_19W
cd /scratch/sf040090/genebycell_19W

echo "Now I will sort and index the transcriptome files"
samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/25transcriptome/25_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 25_txn_featureassigned_sorted.bam
samtools index 25_txn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/26transcriptome/26_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 26_txn_featureassigned_sorted.bam
samtools index 26_txn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/27transcriptome/27_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 27_txn_featureassigned_sorted.bam
samtools index 27_txn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/28transcriptome/28_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 28_txn_featureassigned_sorted.bam
samtools index 28_txn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/29transcriptome/29_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 29_txn_featureassigned_sorted.bam
samtools index 29_txn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_txn_19W/30transcriptome/30_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 30_txn_featureassigned_sorted.bam
samtools index 30_txn_featureassigned_sorted.bam
echo "The transcriptome files are sorted and indexed"

echo "Now I will sort and index the transposonome files"
samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/25transposonome/25_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 25_tspn_featureassigned_sorted.bam
samtools index 25_tspn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/26transposonome/26_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 26_tspn_featureassigned_sorted.bam
samtools index 26_tspn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/27transposonome/27_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 27_tspn_featureassigned_sorted.bam
samtools index 27_tspn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/28transposonome/28_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 28_tspn_featureassigned_sorted.bam
samtools index 28_tspn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/29transposonome/29_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 29_tspn_featureassigned_sorted.bam
samtools index 29_tspn_featureassigned_sorted.bam

samtools sort /home/sf040090/Li-2017-hFGC/cln_tspn_19W/30transposonome/30_tspn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 30_tspn_featureassigned_sorted.bam
samtools index 30_tspn_featureassigned_sorted.bam
echo "The transposonome files are sorted and indexed"

echo "next, gonna create cell by gene matrices for each file"
echo "transposonome matrix files"
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 25_txn_featureassigned_sorted.bam -S 25txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 26_txn_featureassigned_sorted.bam -S 26txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 27_txn_featureassigned_sorted.bam -S 27txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 28_txn_featureassigned_sorted.bam -S 28txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 29_txn_featureassigned_sorted.bam -S 29txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 30_txn_featureassigned_sorted.bam -S 30txn_counts.tsv.gz

echo "transcriptome matrix files"
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 25_tspn_featureassigned_sorted.bam -S 25tspn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 26_tspn_featureassigned_sorted.bam -S 26tspn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 27_tspn_featureassigned_sorted.bam -S 27tspn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 28_tspn_featureassigned_sorted.bam -S 28tspn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 29_tspn_featureassigned_sorted.bam -S 29tspn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 30_tspn_featureassigned_sorted.bam -S 30tspn_counts.tsv.gz
echo "cell by gene matrices for each file are DONE:D"

mv /scratch/sf040090/genebycell_19W /home/sf040090/Li-2017-hFGC