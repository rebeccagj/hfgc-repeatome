#!/bin/sh
## Set the job name
#PBS -N babadook
#PBS -l nodes=1:ppn=4,vmem=16gb
# Run my job

mkdir /scratch/sf040090/clean

# 25
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199325_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199325_2.fastq.gz out1=/scratch/sf040090/clean/25_1_clean.fastq.gz out2=/scratch/sf040090/clean/25_2_clean.fastq.gz qtrim=rl trimq=10

# 26
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199326_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199326_2.fastq.gz out1=/scratch/sf040090/clean/26_1_clean.fastq.gz out2=/scratch/sf040090/clean/26_2_clean.fastq.gz qtrim=rl trimq=10

# 27
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199327_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199327_2.fastq.gz out1=/scratch/sf040090/clean/27_1_clean.fastq.gz out2=/scratch/sf040090/clean/27_2_clean.fastq.gz qtrim=rl trimq=10

# 28
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199328_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199328_2.fastq.gz out1=/scratch/sf040090/clean/28_1_clean.fastq.gz out2=/scratch/sf040090/clean/28_2_clean.fastq.gz qtrim=rl trimq=10

# 29
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199329_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199329_2.fastq.gz out1=/scratch/sf040090/clean/29_1_clean.fastq.gz out2=/scratch/sf040090/clean/29_2_clean.fastq.gz qtrim=rl trimq=10

# 30
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g in1=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199330_1.fastq.gz in2=/home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199330_2.fastq.gz out1=/scratch/sf040090/clean/30_1_clean.fastq.gz out2=/scratch/sf040090/clean/30_2_clean.fastq.gz qtrim=rl trimq=10

mv /scratch/sf040090/clean /home/sf040090

# http://seqanswers.com/forums/showthread.php?t=42776 for reference
# oh I'm dumb, I need to assign the cell ID and UMI before cleaning the reads
# however what a delightful trial run this was