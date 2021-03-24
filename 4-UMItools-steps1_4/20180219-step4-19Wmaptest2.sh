#!/bin/sh
## Set the job name
#PBS -N step4-19W
#PBS -l nodes=1:ppn=8,vmem=60gb
# Run my job

mkdir /scratch/sf040090/19Wmap2

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/19Wextract/25_1_extracted.fastq.gz /home/sf040090/19Wextract/25_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0.3 \
     --outFilterMatchNminOverLread 0.3;

mv /scratch/sf040090/19Wmap2 /home/sf040090

# these updates https://github.com/alexdobin/STAR/issues/169

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5324221/ reference for these STAR options