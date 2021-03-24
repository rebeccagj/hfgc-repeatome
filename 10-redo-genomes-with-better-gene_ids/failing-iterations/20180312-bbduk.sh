#!/bin/sh
## Set the job name
#PBS -N 20180312-bbduk
#PBS -l nodes=1:ppn=4,vmem=16gb
# Run my job

# 25
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/25_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/25_#_clean.fastq.gz \
qtrim=rl trimq=10

# 26
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/26_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/26_#_clean.fastq.gz \
qtrim=rl trimq=10

# 27
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/27_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/27_#_clean.fastq.gz \
qtrim=rl trimq=10

# 28
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/28_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/28_#_clean.fastq.gz \
qtrim=rl trimq=10

# 29
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/29_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/29_#_clean.fastq.gz \
qtrim=rl trimq=10

# 30
/home/sf040090/software/bbmap/bbduk.sh -Xmx15g \
in=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wextract/30_#_extracted.fastq.gz \
out=/home/sf040090/Li-2017-hFGC/fastq/male/19W/19Wclean/30_#_clean.fastq.gz \
qtrim=rl trimq=10

echo "20180312-bbduk done"

# http://seqanswers.com/forums/showthread.php?t=42776

# ugh this is literally the dumbest thing. i can get my job to run when i prototype in the in CL but not as a submitted script. ¯\_(ツ)_/¯
# see error log. ran them manually because i hate myself