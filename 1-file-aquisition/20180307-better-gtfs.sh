# reference: http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format

# With downloaded refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/genePredToGtf
chmod +x genePredToGtf
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
gunzip refGene.txt
cut -f 2- refGene.txt > refGene.input
./genePredToGtf file hg38_UCSC-refGene.input hg38_UCSC-refGene.gtf
cat hg38_UCSC-refGene.gtf | sort -k1,1 -k4,4n > hg38_UCSC_refGene_sorted.gtf

# still doesn't have the chrM mito genes. rolling my eyes. to get that...
# UCSC table browser settings: mammal > human > dec 2013 (GRCg38/hg38) > Genes and Gene Predictions > All GENCODE V27 > Basic (wgEncodeGencodeBasicV27) > position chrM:1-16569 > output: all fields from selected table > output file chrM.txt
cut -f 2- chrM.txt > chrM.input
emacs chrM.input # opens text editor so you can delete the row one column name information
./genePredToGtf file chrM.input chrM.gtf
cat chrM.gtf | sort -k1,1 -k4,4n > chrM_sorted.gtf

# bring transposon the files together
wc -l hg38_UCSC_refGene_sorted.gtf # 1471163
wc -l chrM_sorted.gtf # 107
cat chrM_sorted.gtf hg38_UCSC_refGene_sorted.gtf > hg38_UCSC_refGene_chrM.gtf
wc -l hg38_UCSC_refGene_chrM.gtf # 1471270
cat hg38_UCSC_refGene_chrM.gtf | sort -k1,1 -k4,4n > hg38_UCSC_refGene_chrM_sorted.gtf

#  transposonome 
# https://drive.google.com/uc?export=download&confirm=ZmK_&id=0B4W3iBaKAOm2NDVKbGh3VFlLeWs
grep -v "(*)n" hg38.fa.gtf > hg38_UCSC_repeatmasker_noSR.gtf # removed the simple repeats
cat hg38_UCSC_repeatmasker_noSR.gtf | sort -k1,1 -k4,4n > hg38_UCSC_repeatmasker_noSR_sorted.gtf

# these new files were added the the genomefiles directory and the README.txt updated