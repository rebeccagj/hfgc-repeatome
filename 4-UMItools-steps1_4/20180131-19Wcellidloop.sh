#!/bin/sh
## Set the job name
#PBS -N umiloop
#PBS -l nodes=1:ppn=4,vmem=32gb
# Run my job

FILES=/home/sf040090/Li-2017-hFGC/fastq/male/19W/*_2.fastq.gz
for f in $FILES; do \
    umi_tools whitelist --stdin "$f" --bc-pattern=CCCCCCCCNNNNNNNN --expect-cells=60 --log2stderr > "$f"-whitelist.txt;
done

# generates whitelist.txt files for everything in the loop
#  18 SRR4199325_2.fastq.gz-whitelist.txt
#  24 SRR4199326_2.fastq.gz-whitelist.txt
#  48 SRR4199327_2.fastq.gz-whitelist.txt
#  19 SRR4199328_2.fastq.gz-whitelist.txt
#  48 SRR4199329_2.fastq.gz-whitelist.txt
#  17 SRR4199330_2.fastq.gz-whitelist.txt
# 174 total

# going to rerun with specified number of cells from LiLi et al
# this files are stashed in expectcells60-whitelists