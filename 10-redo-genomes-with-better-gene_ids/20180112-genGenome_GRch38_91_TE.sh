# generating a new reference genome for Li et al STAR alignment

# download the most recent UNMASKED GRCh38 FASTA files
mkdir GRCh38/GRch38-91
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > ensemblGRCh38_r91_unmasked.fa

# initialize CBC star module
module load CBC
module load star

# generate reference genome
STAR --runThreadN 11 --runMode genomeGenerate --genomeDir /home/sf040090/Li-2017-hFGC/hg38_79-TE-150overhang --genomeFastaFiles  /home/sf040090/Li-2017-hFGC/GRCh38_r91_unmasked.fa --sjdbGTFfile /home/sf040090/Li-2017-hFGC/hg38.fa.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 31000000000
# runThread = # threads requested
# genomeDir = a new dir made to store the genome generated
# sjdbOverhang = sjdbOverhang set as readlength -1. So if you have 75 bp read then it should be set to 74
# limitGenomeGenerateRAM 31000000000 = 31 gig limit