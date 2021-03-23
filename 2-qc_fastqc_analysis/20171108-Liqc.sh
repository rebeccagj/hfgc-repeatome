#!/bin/sh
## Set the job name
#PBS -N 20171108_Liqc
#PBS -l nodes=1:ppn=32,vmem=64gb
# Run my job

mkdir /scratch/sf040090/20171108-Li-qc
cd /scratch/sf040090/20171108-Li-qc

module load CBC
module load fastqc

fastqc /home/sf040090/Li-2017-hFGC/fastq/SRR*

mv /scratch/sf040090/20171108-Li-qc /home/sf040090/