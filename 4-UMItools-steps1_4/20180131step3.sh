#!/bin/sh
## Set the job name
#PBS -N 20180131step3
#PBS -l nodes=1:ppn=8,vmem=32gb
# Run my job

mkdir /scratch/sf040090/SRR4199325scratch

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
--stdin /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_2.fastq.gz \
--stdout /scratch/sf040090/SRR4199325scratch/SRR4199325_2_extracted.fastq.gz \
--read2-in /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_1.fastq.gz \
--read2-out /scratch/sf040090/SRR4199325scratch/SRR4199325_1_extracted.fastq.gz \
--filter-cell-barcode \
--whitelist /home/sf040090/SRR4199325_2whitelist.txt

mv /scratch/sf040090/SRR4199325scratch /home/sf040090

# this job tests up through step 3
