#!/bin/sh
## Set the job name
#PBS -N 20180311-ensembl_transcriptome
#PBS -l nodes=1:ppn=4,vmem=60gb
# Run my job

mkdir /scratch/sf040090/transcriptome_CRCh38_91_100overhang
cd /scratch/sf040090/transcriptome_CRCh38_91_100overhang

module load CBC
module load star/2.5.4b 

STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir /scratch/sf040090/transcriptome_CRCh38_91_100overhang \
--genomeFastaFiles /home/sf040090/Li-2017-hFGC/genomefiles/ensemblGRCh38_r91_unmasked.fa \
--sjdbGTFfile /home/sf040090/Li-2017-hFGC/genomefiles/ensemblGRCh38_r91.gtf \
--sjdbOverhang 100 \
--limitGenomeGenerateRAM 25000000000 \
--limitSjdbInsertNsj 370000

mv /scratch/sf040090/transcriptome_CRCh38_91_100overhang /home/sf040090/Li-2017-hFGC/

echo "20180310-ensembl_transcriptome finished"

# ugh finally found this link: https://groups.google.com/forum/#!topic/rna-star/ddhJDgvZfNA
# this problem is caused by the genome+junctions length being close to the 2^32 boundary.
# To avoid this problem you would need to specify more precisely the total number of 
# (collapsed) junctions to be inserted in the genome. 
# You can find this number in the Log.out file:
# Finished SA search: number of new junctions=3292413, old junctions=0
# If you specify --limitSjdbInsertNsj 3400000, the error should go away.

# YAY! THIS WORKED! Mar 12 12:02:07 ..... finished successfully