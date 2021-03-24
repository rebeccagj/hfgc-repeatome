#!/bin/sh
## Set the job name
#PBS -N 20180307-new_gene_id_genomes
#PBS -l nodes=1:ppn=1,vmem=160gb
# Run my job

mkdir /scratch/sf040090/genomes20180307
cd /scratch/sf040090/genomes20180307
mkdir transcriptome_hg38_150overhang
mkdir transposonome_hg38_150overhang

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180307/transcriptome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_refGene_chrM_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 150000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/Log.out /home/sf040090/transcriptomeLog.out

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180307/transposonome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_repeatmasker_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 150000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/Log.out /home/sf040090/transposonomeLog.out
mv /scratch/sf040090/genomes20180307 /home/sf040090/Li-2017-hFGC

echo "genomegen20180307 finished"

# ran on n13

# job hung at the "sorting Suffix Array chunks and saving them to disk" step, per typical
# increased job threads to 4 per this thread (https://github.com/alexdobin/STAR/issues/202) before trying to modify the limitGenomeGenerateRAM coordinate