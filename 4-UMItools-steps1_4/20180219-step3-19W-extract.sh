#!/bin/sh
## Set the job name
#PBS -N step3-19W-28_30
#PBS -l nodes=1:ppn=8,vmem=40gb
# Run my job

mkdir /scratch/sf040090/19Wextract28_30

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199328_2p.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract28_30/28_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199328_1p.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract28_30/28_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/28_2-whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199329_2p.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract28_30/29_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199329_1p.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract28_30/29_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/29_2-whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199330_2p.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract28_30/30_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199330_1p.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract28_30/30_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/30_2-whitelist.txt


mv /scratch/sf040090/19Wextract28_30 /home/sf040090


