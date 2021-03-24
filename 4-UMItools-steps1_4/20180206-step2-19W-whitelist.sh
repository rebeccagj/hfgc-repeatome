#!/bin/sh
## Set the job name
#PBS -N step2-19W
#PBS -l nodes=1:ppn=4,vmem=32gb
# Run my job

mkdir /scratch/sf040090/19Wwhitelist/

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199325_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=24 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/25_2-whitelist.txt;

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199326_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=25 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/26_2-whitelist.txt;

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199327_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=48 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/27_2-whitelist.txt;

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199328_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=25 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/28_2-whitelist.txt;

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199329_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=50 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/29_2-whitelist.txt;

umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199330_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=24 \
                    --log2stderr > /scratch/sf040090/19Wwhitelist/30_2-whitelist.txt;

mv /scratch/sf040090/19Wwhitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W

# ran on n12