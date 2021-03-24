#!/bin/sh
## Set the job name
#PBS -N 20180308a-new_gene_id_genomes
#PBS -l nodes=1:ppn=4,vmem=160gb
# Run my job

mkdir /scratch/sf040090/genomes20180308a
cd /scratch/sf040090/genomes20180308a
mkdir transcriptome_hg38_150overhang
mkdir transposonome_hg38_150overhang

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180308a/transcriptome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_refGene_chrM_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 25000000000 \
--limitSjdbInsertNsj 6000000

# EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem
# SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star
# am going to update STAR as well as --sjdbOverhang 100 per this thread https://groups.google.com/forum/#!topic/rna-star/PCX4IviXch0

cp /scratch/sf040090/genomes20180308a/Log.out /home/sf040090/transcriptomeLog.out

/home/sf040090/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/genomes20180308a/transposonome_hg38_150overhang \
--genomeFastaFiles  /home/sf040090/Li-2017-hFGC/genomefiles/UCSChg38.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/hg38_UCSC_repeatmasker_sorted.gtf \
--sjdbOverhang 149 \
--limitGenomeGenerateRAM 250000000000 \
--limitSjdbInsertNsj 6000000

cp /scratch/sf040090/genomes20180308a/Log.out /home/sf040090/transposonomeLog.out
mv /scratch/sf040090/genomes20180308a /home/sf040090/Li-2017-hFGC

echo "genomegen20180308 finished"

# same as "20180308-new_gene_id_genomes.sh" except with --limitGenomeGenerateRAM 250000000000 per this thread: https://github.com/alexdobin/STAR/issues/364
# did not complete the transcriptome, although did finish the transposonome.
# transcriptome error still: EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem
SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star
# per this thread, https://groups.google.com/forum/#!topic/rna-star/PCX4IviXch0 i updated sjdboverhang and (and irrelevantly also increased my runThread)