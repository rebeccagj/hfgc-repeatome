#!/bin/sh
## Set the job name
#PBS -N 19Wmap26_30
#PBS -l nodes=1:ppn=8,vmem=60gb
# Run my job

mkdir /scratch/sf040090/19Wmap26_30
cd /scratch/sf040090/19Wmap26_30

mkdir 26
cd 26
/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/26_1_extracted.fastq.gz /home/sf040090/19Wextract/26_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap26_30

mkdir 27
cd 27
/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/27_1_extracted.fastq.gz /home/sf040090/19Wextract/27_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap26_30

mkdir 28
cd 28
/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/28_1_extracted.fastq.gz /home/sf040090/19Wextract/28_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap26_30

mkdir 29
cd 29
/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/29_1_extracted.fastq.gz /home/sf040090/19Wextract/29_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap26_30

mkdir 30
cd 30
/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/30_1_extracted.fastq.gz /home/sf040090/19Wextract/30_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap26_30


mv /scratch/sf040090/19Wmap26_30 /home/sf040090

# update = finishing the alignment of the other 5 19W fastq sets