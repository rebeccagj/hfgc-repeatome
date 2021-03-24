# better gtf analysis! (:)
setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/genebycell_19W/") #set wd to directory with UMI_Tools tsv file outputs

#function for inputting data, appending column name with SRR ID number, and assigning to global environment ####
library(data.table)
readrenamecol = function(arg1, arg2, arg3){
  arg1 = read.delim(file = arg2, header = T, sep = "\t", row.names = 1)
  colnames(arg1)[1:ncol(arg1)] = paste(colnames(arg1)[1:ncol(arg1)], arg3, sep = "_")
  arg1 = data.table(arg1, keep.rownames = TRUE) #data.frame objects are too large, using data.table object instead
  return(arg1)
}
#don't worry this is going to have genes as column one, rather than row one.
# you'll take care of that when you read the data in

cmncol = read.table("both/old/commongenes.txt", row.names = 1)

# cat te/txn combo files ####
#read in combined tspn and txn datasets for W19  
tspntxn25 = readrenamecol(tspntxn25, "both/25_tspntxn.tsv", "25")
tspntxn26 = readrenamecol(tspntxn26, "both/26_tspntxn.tsv", "26")
tspntxn27 = readrenamecol(tspntxn27, "both/27_tspntxn.tsv", "27")
tspntxn28 = readrenamecol(tspntxn28, "both/28_tspntxn.tsv", "28")
tspntxn29 = readrenamecol(tspntxn29, "both/29_tspntxn.tsv", "29")
tspntxn30 = readrenamecol(tspntxn30, "both/30_tspntxn.tsv", "30")
#it's ok future Rebecca, different numbers of obs. comes from discrete cell-to-cell, fastq-to-fastq alignments

#concatinate files, dropping rows that are not present in all SRR sets
txn_te_19W = merge(counts25, counts26)
txn_te_19W = merge.data.frame(txn_te_19W, counts27)
txn_te_19W = merge.data.frame(txn_te_19W, counts28)
txn_te_19W = merge.data.frame(txn_te_19W, counts29)
txn_te_19W = merge.data.frame(txn_te_19W, counts30)
write.table(txn_te_19W, file = "txn_te_19W.txt", quote = FALSE, sep = "\t", row.names = FALSE) #write merged data.table to txt file
#HAHAHHA IT'S FUCKING DONE
#also i changed 1:1 from "rn" to "genes" in the txt file manually via a text editor

# cat txn alone files ####
#read in txn alone datasets for W19
txn25 = readrenamecol(txn25, "txn/25txn_counts.tsv", "25")
txn26 = readrenamecol(txn26, "txn/26txn_counts.tsv", "26")
txn27 = readrenamecol(txn27, "txn/27txn_counts.tsv", "27")
txn28 = readrenamecol(txn28, "txn/28txn_counts.tsv", "28")
txn29 = readrenamecol(txn29, "txn/29txn_counts.tsv", "29")
txn30 = readrenamecol(txn30, "txn/30txn_counts.tsv", "30")
#it's ok future Rebecca, different numbers of obs. comes from discrete cell-to-cell, fastq-to-fastq alignments

#concatinate all datasets, dropping rows that are not present in all SRR sets
txn_19W = merge(txn25, txn26)
txn_19W = merge(txn_19W, txn27)
txn_19W = merge(txn_19W, txn28)
txn_19W = merge(txn_19W, txn29)
txn_19W = merge(txn_19W, txn30)
write.table(txn_19W, file = "txn_19W.txt", quote = FALSE, sep = "\t", row.names = FALSE) #write merged data.table to txt file
#also i changed 1:1 from "rn" to "genes" in the txt file manually via a text editor

# cat te alone files before this code stops working ####
#read in tspn alone datasets for W19
te25 = readrenamecol(te25, "tspn/25tspn_counts_nooverlap.tsv", "25")
te26 = readrenamecol(te26, "tspn/26tspn_counts_nooverlap.tsv", "26")
te27 = readrenamecol(te27, "tspn/27tspn_counts_nooverlap.tsv", "27")
te28 = readrenamecol(te28, "tspn/28tspn_counts_nooverlap.tsv", "28")
te29 = readrenamecol(te29, "tspn/29tspn_counts_nooverlap.tsv", "29")
te30 = readrenamecol(te30, "tspn/30tspn_counts_nooverlap.tsv", "30")
#it's ok future Rebecca, different numbers of obs. comes from discrete cell-to-cell, fastq-to-fastq alignments

#concatinate all datasets, dropping rows that are not present in all SRR sets
te_19W = merge(te25, te26)
te_19W = merge(te_19W, te27)
te_19W = merge(te_19W, te28)
te_19W = merge(te_19W, te29)
te_19W = merge(te_19W, te30)
write.table(te_19W, file = "te_19W.txt", quote = FALSE, sep = "\t", row.names = FALSE) #write merged data.table to txt file
#also i changed 1:1 from "rn" to "genes" in the txt file manually via a text editor

#whole week 19 datasets ####
txn_te_19W = read.delim(file = "txn_te_19W.txt", row.names = 1) 
txn_19W = read.delim(file = "txn_19W.txt", row.names = 1)
te_19W = read.delim(file = "te_19W.txt", row.names = 1)

setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/GSE86146_RAW/") #set wd to directory with UMI_Tools tsv file outputs
pub19W = read.csv("GSE86146_RAW_all19Wcounts.txt", sep = "\t", row.names = 2)[-1]

setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/genebycell_19W/") #set wd to directory with UMI_Tools tsv file outputs
library(Seurat)
txn19W_matrix = as(as.matrix(txn_19W), "dgCMatrix")
w19 <- CreateSeuratObject(raw.data = txn19W_matrix, min.cells = 5, is.expr = 1, project = "hFGC19")

mitogenes <- grep(pattern = "^MT-", x = rownames(x = w19@data), value = TRUE)
percentmito <- Matrix::colSums(w19@raw.data[mitogenes, ])/Matrix::colSums(w19@raw.data)
w19 <- AddMetaData(object = w19, metadata = percentmito, col.name = "percentmito")
length(x = w19@cell.names) #196

par(mfrow = c(1, 2))
GenePlot(object = w19, gene1 = "nUMI", gene2 = "percentmito")
GenePlot(object = w19, gene1 = "nUMI", gene2 = "nGene")
#w19 = FilterCells(object = w19, subset.names = c("nGene", "percentmito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.07))
#decided not to filter

w19 <- NormalizeData(object = w19, normalization.method = "LogNormalize", scale.factor = 10000)
w19 <- FindVariableGenes(object = w19, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = Inf, y.cutoff = 0.5)
length(x = w19@var.genes) #804

w19 <- ScaleData(object = w19, vars.to.regress = c("nUMI", "percentmito"))
w19 <- RunPCA(object = w19, pc.genes = w19@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = w19, pcs.use = 1:5)
PCAPlot(object = w19, 1,2)
PCAPlot(object = w19, 1,3)
PCAPlot(object = w19, 3,2)
PCHeatmap(object = w19, pc.use = 1, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = w19, pc.use = 2, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = w19, pc.use = 3, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
w19 <- JackStraw(object = w19, num.replicate = 50, do.print = FALSE)
JackStrawPlot(object = w19, PCs = 1:12) #1:9 significant
PCElbowPlot(object = w19)
w19 <- FindClusters(object = w19, reduction.type = "pca", dims.use = 1:9, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = w19)
w19 <- RunTSNE(object = w19, dims.use = 1:9, do.fast = TRUE)
TSNEPlot(object = w19)
save(w19, file = "~/20180321-seuratw19.Robj")
