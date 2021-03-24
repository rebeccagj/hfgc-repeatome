#!/bin/sh
## Set the job name
#PBS -N 20180308-new_gene_id_genomes
#PBS -l nodes=1:ppn=4,vmem=160gb
# Run my job

mkdir /scratch/sf040090/genomes20180308
cd /scratch/sf040090/genomes20180308
mkdir transcriptome_hg38_150overhang
mkdir transposonome_hg38_150overhang

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180308/transcriptome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_refGene_chrM_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 150000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/genomes20180308/Log.out /home/sf040090/transcriptomeLog.out

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180308/transposonome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_repeatmasker_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 150000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/genomes20180308/Log.out /home/sf040090/transposonomeLog.out
mv /scratch/sf040090/genomes20180308 /home/sf040090/Li-2017-hFGC

echo "genomegen20180308 finished"

# job 1083832 failed because i had a path messed up from yesterday's job
# fixed and rerun as job 1083839
# job did not work. modified limitGenomeGenerateRAM and resubmitted as 20180308a