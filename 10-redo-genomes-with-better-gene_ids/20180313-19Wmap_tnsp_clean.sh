#!/bin/sh
## Set the job name
#PBS -N 20180313-19Wmap_tnsp_clean
#PBS -l nodes=1:ppn=6,vmem=60gb
# Run my job

module load CBC
module load star/2.5.3a

mkdir /scratch/sf040090/19Wmap20180313
cd /scratch/sf040090/19Wmap20180313

mkdir 25transposonome
cd 25transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/25_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/25_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 26transposonome
cd 26transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/26_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/26_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 27transposonome
cd 27transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/27_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/27_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 28transposonome
cd 28transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/28_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/28_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 29transposonome
cd 29transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/29_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/29_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 30transposonome
cd 30transposonome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transposonome_hg38_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/30_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/30_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

echo "20180313-19Wmap_tnsp_clean complete"

# redoing alignment of W19 with better gtfs and cleaned reads
# successful run!