#!/bin/sh
## Set the job name
#PBS -N step4-19W
#PBS -l nodes=1:ppn=8,vmem=32gb
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

# EXITING: fatal error trying to allocate genome arrays, exception thrown: std::bad_alloc
# Possible cause 1: not enough RAM. Check if you have enough RAM 49496831129 bytes
# Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v 49496831129
# running with 60gb ram requested next


