#!/bin/sh
## Set the job name
#PBS -N 28_30pipe
#PBS -l nodes=1:ppn=1,vmem=32gb
# Run my job

mkdir /scratch/sf040090/pipe-edits

gunzip -c /home/sf040090/Li-2017-hFGC/fastq/male/19W/stupidfiles/SRR4199328_1.fastq.gz > /scratch/sf040090/pipe-edits/SRR4199328_1.fastq
gunzip -c /home/sf040090/Li-2017-hFGC/fastq/male/19W/stupidfiles/SRR4199328_2.fastq.gz > /scratch/sf040090/pipe-edits/SRR4199328_2.fastq
gunzip -c /home/sf040090/Li-2017-hFGC/fastq/male/19W/stupidfiles/SRR4199329_1.fastq.gz > /scratch/sf040090/pipe-edits/SRR4199329_1.fastq
gunzip -c /home/sf040090/Li-2017-hFGC/fastq/male/19W/stupidfiles/SRR4199329_2.fastq.gz > /scratch/sf040090/pipe-edits/SRR4199329_2.fastq

cd /scratch/sf040090/pipe-edits

cat SRR4199328_1.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199328_1p.fastq
cat SRR4199328_2.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199328_2p.fastq
cat SRR4199329_1.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199329_1p.fastq
cat SRR4199329_2.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199329_2p.fastq

gzip /scratch/sf040090/pipe-edits/SRR41993*

mv /scratch/sf040090/pipe-edits/ /home/sf040090/Li-2017-hFGC/fastq/male/19W/