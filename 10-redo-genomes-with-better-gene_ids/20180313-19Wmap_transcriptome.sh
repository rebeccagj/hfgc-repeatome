#!/bin/sh
## Set the job name
#PBS -N 20180313-19Wmap_transcriptome
#PBS -l nodes=1:ppn=6,vmem=80gb
# Run my job

module load CBC
module load star/2.5.3a 

mkdir /scratch/sf040090/19Wmap20180313
cd /scratch/sf040090/19Wmap20180313

mkdir 25transcriptome
cd 25transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/25_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/25_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 26transcriptome
cd 26transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/26_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/26_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 27transcriptome
cd 27transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/27_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/27_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 28transcriptome
cd 28transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/28_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/28_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 29transcriptome
cd 29transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/29_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/29_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

mkdir 30transcriptome
cd 30transcriptome
STAR --runThreadN 6 \
     --genomeDir /home/sf040090/Li-2017-hFGC/transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/30_1_clean.fastq.gz /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/30_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
cd /scratch/sf040090/19Wmap20180313

echo "20180313-19Wmap_transcriptome complete"

# successful run