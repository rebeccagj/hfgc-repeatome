#I got this code from here: https://www.biostars.org/p/93494/
#using it to download all the SRA files located here: https://www.ncbi.nlm.nih.gov/sra?term=SRP083134
#or here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP083134&go=go
#also this reddit thread was v helpful https://tinyurl.com/ybd3x983

#load bioconductor
source('http://bioconductor.org/biocLite.R')
#get the package we want
biocLite('SRAdb')
bioLite('DBI')
a
y

#load the package
library(SRAdb)
library(DBI)

#shorthand
srafile = getSRAdbFile() #this may take a long time ughhhh
con = dbConnect(RSQLite::SQLite(), srafile)

#query the local SQLite database
listSRAfile('SRP083134',con)
#there are 43 SRAs to download

#determine the location for downloading the SRAs
setwd("/Volumes/DATA$/Users/Rebecca-Jaszczak/public/single-cell_datasets/Li-etal-2017/SRA")

#have R do the downloads for you
#srcType can be ftp or fasp; fasp is better for larger files
#ascpCMD: download the free ascp connect client http://downloads.asperasoft.com/en/downloads/8?list
#follow the ascpCMD setup instructions in the help file
help(getSRAfile)
getSRAfile('SRP083134', con, fileType='sra', srcType = 'fasp', ascpCMD = "'/Users/Rebecca/Applications/Aspera-Connect.app/Contents/Resources/ascp' -QT -l 300m -i '/Users/Rebecca/Applications/Aspera-Connect.app/Contents/Resources/asperaweb_id_dsa.putty'")
