[toc]

# Analysis Pipeline for 19W Supplemental Data

assumes the following software packages are installed: umi_tools (0.5.4), STAR(2.5.4b), BBTools, samtools(1.6), Subread (1.6.0)

(12/2018 note) I actually recommend updating samtools to the most recent version, because they implimented options to throw more threads and more RAM at your command, which is really helpful for these large files.

##  1. FIND CELL BARCODE WHITELIST
```
umi_tools whitelist --stdin SRR4199325_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=24 \
                    --log2stderr > 25_2-whitelist.txt
umi_tools whitelist --stdin SRR4199326_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=25 \
                    --log2stderr > 26_2-whitelist.txt
umi_tools whitelist --stdin SRR4199327_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=48 \
                    --log2stderr > 27_2-whitelist.txt
umi_tools whitelist --stdin SRR4199328_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=25 \
                    --log2stderr > 28_2-whitelist.txt
umi_tools whitelist --stdin SRR4199329_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=50 \
                    --log2stderr > 29_2-whitelist.txt
umi_tools whitelist --stdin SRR4199330_2.fastq.gz \
                    --bc-pattern=CCCCCCCCNNNNNNNN \
                    --set-cell-number=24 \
                    --log2stderr > 30_2-whitelist.txt
```

I ran into an issue. The pipelines accepts read SRR4199325.1 1/ **1** and SRR4199325.1 1/ __2__ where the bold digit after the / is the fwd/rev pass read. This nomenclature doesn't cause hiccups down the line. However, some SRRs were formatted differently. The pipelines denies formatting SRR4199330.1. **1** 1 and SRR4199330.1. **2** 1 as the ordering of the 1 and 2 fwd/rev pass nomenclature caused mistakes. SRR4199328, SRR4199329, and SRR4199330 had this issue. End characters of these reads were trimmed to output SRR4199330.1 and SRR4199330.1 in the fwd and rev read files, respectively, so that the F and R pairs would match up correctly.

```
gunzip -c SRR4199328_1.fastq.gz > SRR4199328_1.fastq # unzips the files so they can be modifed
gunzip -c SRR4199328_2.fastq.gz > SRR4199328_2.fastq
gunzip -c SRR4199329_1.fastq.gz > SRR4199329_1.fastq
gunzip -c SRR4199329_2.fastq.gz > SRR4199329_2.fastq
gunzip -c SRR4199330_1.fastq.gz > SRR4199330_1.fastq
gunzip -c SRR4199330_2.fastq.gz > SRR4199330_2.fastq

cat SRR4199328_1.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199328_1p.fastq # trims appropriate characters from the read name starting with line 1 and every 5th line following
cat SRR4199328_2.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199328_2p.fastq
cat SRR4199329_1.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199329_1p.fastq
cat SRR4199329_2.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199329_2p.fastq
cat SRR4199330_1.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199330_1p.fastq
cat SRR4199330_2.fastq | awk '{print $1}' | sed '1~2s/..$//g' > SRR4199330_2p.fastq
gzip SRR41993* # recompresses the files for ease of use
```

### From here down, these commands were batched on files 25,26,27,28,29,30, but I am only showing the command run on set 25 as an example.

## 2. EXTRACT CB/UMIS AND FILTER CBS

This command extracts the UMIs from the read they're contained in, and appends them to the name of the read.

```
umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin SRR4199325_2.fastq.gz \
                  --stdout 25_2_extracted.fastq.gz \
                  --read2-in SRR4199325_1.fastq.gz \
                  --read2-out 25_1_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist 25_2-whitelist.txt
```
 Next, using BBTools bbduk, I quality-trim extracted reads to Q10 using Phred algorithm; qtrim=rl trims from left & right.
`bbduk.sh -Xmx15g in=25_#_extracted.fastq.gz out=25_#_clean.fastq.gz qtrim=rl trimq=10` 

## 3. Acquire alignment files and generate STAR reference genomes
In order to align the UMI-assigned, cleaned reads, I acquired the following genome files...

For the transcriptome:

```
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > ensemblGRCh38_r91_unmasked.fa
wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz
gunzip Homo_sapiens.GRCh38.91.gtf.gz
mv Homo_sapiens.GRCh38.91.gtf ensemblGRCh38_r91.gtf
rm Homo_sapiens.GRCh38.dna.chromosome.*
```

For the repeatome:

```
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa hg38_UCSC.fa
# hg38_UCSC_repeatmasker_sorted.gtf - repeatome annotations from repeatmasker hg38.fa.out, converted to a GTF by Ping Ye
grep -v "(*)n" hg38_UCSC_repeatmasker_sorted.gtf > hg38_UCSC_repeatmasker_noSR.gtf # removed the simple repeats
cat hg38_UCSC_repeatmasker_noSR.gtf | sort -k1,1 -k4,4n > hg38_UCSC_repeatmasker_noSR_sorted.gtf 
```

To align with `star/2.5.4b`, I did a but of research about making the STAR genomes. The `sjdbOverhang` flag should be *readlength minus 1*, but Alex Dobin (STAR creator) says >100 doesn't improve accuracy. https://groups.google.com/d/msg/rna-star/h9oh10UlvhI/ZwAs1xFVCAAJ Alex also recommended I specific `limitGenomeGenerateRAM` as without specifying this variable, jobs fail. Finally, I add `limitSjdbInsertNsj` to more precisely specify the number of collapsed junctions to be inserted in the genome; this is appended on a case by case basis based on failed runs indication I needed to specify more numbers.

```
# repeatome
mkdir repeatome_hg38_100overhang
STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir repeatome_hg38_100overhang \
--genomeFastaFiles  UCSChg38.fa \ 
--sjdbGTFfile hg38_UCSC_repeatmasker_sorted.gtf \
--sjdbOverhang 100 \ 
--limitGenomeGenerateRAM 250000000000 \ 
--limitSjdbInsertNsj 6000000 

# transcriptome
mkdir transcriptome_CRCh38_91_100overhang
STAR --runThreadN 3 \
--runMode genomeGenerate \
--genomeDir transcriptome_CRCh38_91_100overhang \
--genomeFastaFiles ensemblGRCh38_r91_unmasked.fa \
--sjdbGTFfile ensemblGRCh38_r91.gtf \
--sjdbOverhang 100 \
--limitSjdbInsertNsj 3400000
```

## 4. MAP READS

### Map to transcriptome = txn

```
STAR --runThreadN 6 \
     --genomeDir transcriptome_CRCh38_91_100overhang/ \
     --readFilesIn 25_1_clean.fastq.gz 25_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
```
`outFilterScoreMinOverLread` and `outFilterMatchNminOverLread` were added because the fastq files were unsorted; rather than an arduous sort of these 30gb (when zipped) files, I tweaked these parameters to allow reads to match up with their pairs

Renaming and sorting:

`mv Aligned.sortedByCoord.out.bam 25_txn_aligned.sortedbycoord.out.bam` Renames the output bam to distinguish it from the other 12 bams.
`samtools view 25_txn_aligned.sortedbycoord.out.bam | cut -f 2 | sort | uniq -c` Lists the unique sam flags contained in the alignment
`samtools view -h -f 0x2 -b 25_txn_aligned.sortedbycoord.out.bam > 25_txn_aligned.sortedbycoord.f2.bam` Keeps the header and the reads that are aligned in proper pairs
`samtools view 25_txn_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c` Verifies that the flags maintained are ones which are associated with properly paired reads


### Map to repeatome = rptm

```
STAR --runThreadN 6 \
     --genomeDir repeatome_hg38_100overhang/ \
     --readFilesIn 25_1_clean.fastq.gz 25_2_clean.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 100 \
     --outSAMmultNmax 100 \
     --outFilterMismatchNmax 3 \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0;
```
Flags `outFilterMultimapNmax 100`, `winAnchorMultimapNmax 100`, `outSAMmultNmax 100`, and `outFilterMismatchNmax 3` were modified to maintain multimappers.

Renaming and sorting:

`mv Aligned.sortedByCoord.out.bam 25_rptm_aligned.sortedByCoord.out.bam`
`samtools view 25_rptm_aligned.sortedByCoord.out.bam | cut -f 2 | sort | uniq -c`
`samtools view -hf 2 25_rptm_aligned.sortedByCoord.out.bam > 25_rptm_aligned.sortedbycoord.f2.bam` Note that multi-mappers are kept
`samtools view 25_rptm_aligned.sortedbycoord.f2.bam | cut -f 2 | sort | uniq -c` Verifies that the flags maintained are ones which are associated with properly paired reads, single or multimappers

## 5. ASSIGN READS TO GENES

`featureCounts -a ensemblGRCh38_r91.gtf -o 25_genes_assigned -R BAM 25_txn_aligned.sortedbycoord.f2.bam -T 8 -p -g gene_name` Output is `25_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam`
`samtools sort 25_txn_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 25_txn_featureassigned_sorted.bam`

`featureCounts -a hg38_UCSC_repeatmasker_noSR_sorted.gtf -o 25_rptm_assigned -R BAM 25_rptm_aligned.sortedbycoord.f2.bamM -T 8 -p -g gene_id -MO` where MO is added to continue maintaining multimappers
`samtools sort 25_rptm_aligned.sortedbycoord.f2.bam.featureCounts.bam -o 25_rptm_featureassigned_sorted.bam`

## 6. COUNT UNIQUE READS PER GENES PER CELL

```
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 25_txn_aligned.sortedbycoord.f2.bam -S 25txn_counts.tsv.gz
umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I 25_rptm_featureassigned_sorted.bam -S 25rptm_counts.tsv.gz
```
