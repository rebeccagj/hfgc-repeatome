#!/bin/sh
## Set the job name
#PBS -N step3-19W
#PBS -l nodes=1:ppn=8,vmem=60gb
# Run my job

mkdir /scratch/sf040090/19Wextract

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199325_2.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract/25_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199325_1.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract/25_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/25_2-whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199326_2.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract/26_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199326_1.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract/26_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/26_2-whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199327_2.fastq.gz \
                  --stdout /scratch/sf040090/19Wextract/27_2_extracted.fastq.gz \
                  --read2-in /home/sf040090/Li-2017-hFGC/fastq/male/19W/SRR4199327_1.fastq.gz \
                  --read2-out /scratch/sf040090/19Wextract/27_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist /home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wwhitelist/27_2-whitelist.txt

# for some reason 28, 29, and 30 did not work

mv /scratch/sf040090/19Wextract /home/sf040090


