setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/UMIcountwide/") #set wd to directory with UMI_Tools tsv file outputs

readrenamecol = function(arg1, arg2, arg3){
  arg1 = read.delim(file = arg2, header = T, sep = "\t")
  colnames(arg1)[2:ncol(arg1)] = paste(colnames(arg1)[2:ncol(arg1)], arg3, sep = "_")
  return(arg1)
} #function for inputting data, appending column name with SRR ID number, and assigning to global environment

#read in the 6, 19W samples
counts25 = as.data.frame(readrenamecol(counts25, "25_counts.tsv", "25"))
counts26 = as.data.frame(readrenamecol(counts26, "26_counts.tsv", "26"))
counts27 = as.data.frame(readrenamecol(counts27, "27_counts.tsv", "27"))
counts28 = as.data.frame(readrenamecol(counts28, "28_counts.tsv", "28"))
counts29 = as.data.frame(readrenamecol(counts29, "29_counts.tsv", "29"))
counts30 = as.data.frame(readrenamecol(counts30, "30_counts.tsv", "30"))

#concatinate all datasets, dropping rows that are not present in all SRR sets
all19W = merge.data.frame(counts25,counts26)
all19W = merge.data.frame(all19W, counts27)
all19W = merge.data.frame(all19W, counts28)
all19W = merge.data.frame(all19W, counts29)
all19W = merge.data.frame(all19W, counts30)
write.csv(all19W, file = "all19Wcounts.csv") #write merged data.frame to csv file