#!/bin/sh
## Set the job name
#PBS -N 20170929_sra-fhgc
#PBS -l nodes=1:ppn=32,vmem=96gb
# Run my job

mkdir /scratch/sf040090/Li-2017-hFGC-fastq
cd /scratch/sf040090/Li-2017-hFGC-fastq

module load CBC
module load sratoolkit

for SAMPLE in $(cat /home/sf040090/20170926-SRR_Acc_List.txt); 
do
    fastq-dump --split-files --gzip -I -O /scratch/sf040090/Li-2017-hFGC-fastq $SAMPLE
done

cp /scratch/sf040090/Li-2017-hFGC-fastq