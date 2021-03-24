#!/bin/sh
## Set the job name
#PBS -N step4-19W
#PBS -l nodes=1:ppn=8,vmem=60gb
# Run my job

mkdir /scratch/sf040090/19Wmap

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 \
     --genomeDir /home/sf040090/Li-2017-hFGC/UCSC-hg38_rpandtxs-150overhang \
     --readFilesIn /home/sf040090/25_1_extracted.fastq.gz /home/sf040090/25_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate;

mv /scratch/sf040090/19Wextract /home/sf040090

# ran on n13

# RAM allocation worked, and job aligned in just two hours
# however, % of reads unmapped: too short |	91.63%
# 91% unmapped reads in terrible. I think this might be because my fastqs aren't sorted
# https://github.com/alexdobin/STAR/issues/222
# going to sort 25_1 and 25_2 and try again! (:
# https://www.biostars.org/p/15011/

# http://bridgecrest.blogspot.com/2011/08/sort-fastq-file-by-sequence.html
# https://www.biostars.org/p/15011/
