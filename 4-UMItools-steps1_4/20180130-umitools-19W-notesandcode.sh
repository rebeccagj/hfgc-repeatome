# interested in 19W males, based on IF and small RNA seq data
# M_19W_embryo1_101	SRR4199325
# M_19W_embryo1_24	SRR4199326
# M_19W_embryo1_26	SRR4199327
# M_19W_embryo2_102	SRR4199328
# M_19W_embryo2_103	SRR4199329
# M_19W_embryo2_104	SRR4199330

# looking at the file; read 2 contains the barcode info!
# from the supp data of Li et al 2017, we know that there are 8 Cell ID and 8 UMI
zcat SRR4199325_2.fastq.gz | head -n4

# id the real cells
umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_2.fastq.gz --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=25 --log2stderr > SRR4199325_2whitelist.txt

# ugh I got a weird error:
# "ImportError: libgfortran.so.1: cannot open shared object file: 
# No such file or directory" error. upon googling, ran the next lines of code...

conda install libgfortran==1

# now the following code works
umi_tools whitelist --stdin /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_2.fastq.gz --bc-pattern=CCCCCCCCNNNNNNNN --set-cell-number=25 --log2stderr > SRR4199325_2whitelist.txt
# sup table 1 gives me the info on set cell number for M_19W_101 (Li et al found 24 cells)
# this 5G file took ~27mins

wc -l SRR4199325_2whitelist.txt 
# tests length to assure me that there are 25 cells

mkdir /scratch/sf040090/SRR4199325scratch

umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
--stdin /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_2.fastq.gz \
--stdout /scratch/sf040090/SRR4199325scratch/SRR4199325_2_extracted.fastq.gz \
--read2-in /home/sf040090/Li-2017-hFGC/fastq/male/SRR4199325_1.fastq.gz \
--read2-out /scratch/sf040090/SRR4199325scratch/SRR4199325_1_extracted.fastq.gz \
--filter-cell-barcode \
--whitelist /home/sf040090/SRR4199325_2whitelist.txt