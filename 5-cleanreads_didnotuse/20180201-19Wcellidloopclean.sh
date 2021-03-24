#!/bin/sh
## Set the job name
#PBS -N umiloopclean
#PBS -l nodes=1:ppn=4,vmem=32gb
# Run my job

FILES=/home/sf040090/Li-2017-hFGC/fastq/male/19W/clean/*_2_clean.fastq.gz
for f in $FILES; do \
    umi_tools whitelist --stdin "$f" --bc-pattern=CCCCCCCCNNNNNNNN --expect-cells=60 --log2stderr > "$f"-whitelist.txt;
done

# https://github.com/CGATOxford/UMI-tools/issues/213
# this will not run due to cleaned cells not having the bc-pattern specified