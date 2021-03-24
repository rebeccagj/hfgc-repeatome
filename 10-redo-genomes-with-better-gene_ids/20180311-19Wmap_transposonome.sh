#!/bin/sh
## Set the job name
#PBS -N 20180311-19Wmap
#PBS -l nodes=1:ppn=6,vmem=80gb
# Run my job

module load CBC
module load star/2.5.3a 

mkdir /scratch/sf040090/19Wmap20180311
cd /scratch/sf040090/19Wmap20180311

mkdir 25transposonome
cd 25transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/25_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/25_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

mkdir 26transposonome
cd 26transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/26_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/26_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

mkdir 27transposonome
cd 27transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/27_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/27_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

mkdir 28transposonome
cd 28transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/28_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/28_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

mkdir 29transposonome
cd 29transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/29_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/29_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

mkdir 30transposonome
cd 30transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/30_1_extracted.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/30_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180311

echo "20180311-19Wmap complete"

# this worked