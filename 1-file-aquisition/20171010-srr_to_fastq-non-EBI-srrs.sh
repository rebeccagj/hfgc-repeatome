#!/bin/sh
## Set the job name
#PBS -N 20171010-srr_to_fastq-non-EBI-srrs
#PBS -l nodes=1:ppn=18,vmem=60gb
#PBS -M rgjzak@gmail.com
# Run my job

mkdir /scratch/sf040090/srrs-to-fastq
cd /scratch/sf040090/srrs-to-fastq

module load CBC
module load sratoolkit

for SAMPLE in $(cat /home/sf040090/20171010-SRR_Acc_List-SRRs-not-on-EBI.txt); 
do
    fastq-dump --split-files --gzip -I -O /scratch/sf040090/srrs-to-fastq $SAMPLE
done

mv /scratch/sf040090/srrs-to-fastq /home/sf040090/