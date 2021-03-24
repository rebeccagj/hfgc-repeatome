#!/bin/sh
## Set the job name
#PBS -N 20180309-new_gene_id_genomes
#PBS -l nodes=1:ppn=4,vmem=160gb
# Run my job

mkdir /scratch/sf040090/genomes20180309
cd /scratch/sf040090/genomes20180309
mkdir transcriptome_hg38_100overhang
mkdir transposonome_hg38_100overhang

module load CBC
module load star/2.5.4b 

STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180309/transcriptome_hg38_100overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_refGene_chrM_sorted.gtf \
--sjdbOverhang 100 \
--limitGenomeGenerateRAM 25000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/genomes20180309/Log.out /home/sf040090/transcriptomeLog.out

STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180309/transposonome_hg38_100overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_repeatmasker_sorted.gtf \
--sjdbOverhang 100 \
--limitGenomeGenerateRAM 250000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/genomes20180309/Log.out /home/sf040090/transposonomeLog.out
mv /scratch/sf040090/genomes20180309 /home/sf040090/Li-2017-hFGC

echo "genomegen20180308 finished"

# https://groups.google.com/forum/#!topic/rna-star/PCX4IviXch0: updated sjdbOverhang, runThread
# re, Alex on the overlapping or touching exons error: "these messages are just warnings, if you think the GTF file is mostly OK, then STAR will simply ignore the "touching" exons" https://groups.google.com/forum/#!topic/rna-star/pfmAzkEbOm8

# this job succeeds for the transposonome, fails for the transcriptome