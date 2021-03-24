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

#week 19 transcriptome dataset input and clustering ####
txn_te_19W = read.delim(file = "txn_te_19W.txt", row.names = 1) 
txn_19W = read.delim(file = "txn_19W.txt", row.names = 1)
te_19W = read.delim(file = "te_19W.txt", row.names = 1)

setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/GSE86146_RAW/") #set wd to directory with UMI_Tools tsv file outputs
pub19W = read.csv("GSE86146_RAW_all19Wcounts.txt", sep = "\t", row.names = 2)[-1]

setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/genebycell_19W/") #set wd to directory with UMI_Tools tsv file outputs
library(Seurat)
txn19W_matrix = as(as.matrix(txn_19W), "dgCMatrix")

#qc filtering
w19 <- CreateSeuratObject(raw.data = txn19W_matrix, project = "hFGC19")
mitogenes <- grep(pattern = "^MT-", x = rownames(x = w19@data), value = TRUE)
percentmito <- Matrix::colSums(w19@raw.data[mitogenes, ])/Matrix::colSums(w19@raw.data)
w19 <- AddMetaData(object = w19, metadata = percentmito, col.name = "percentmito")
length(x = w19@cell.names) #196
levels(w19@ident) = c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0")
par(mfrow = c(1, 2))
GenePlot(object = w19, gene1 = "nUMI", gene2 = "percentmito", col.use = c("red", "magenta", "darkgreen", "cyan"))
GenePlot(object = w19, gene1 = "nUMI", gene2 = "nGene")
w19 = FilterCells(object = w19, subset.names = c("nGene", "percentmito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(Inf, 0.06))
w19 #An object of class seurat in project hFGC19. 24428 genes across 187 samples.

#normalization and variable genes
w19 <- NormalizeData(object = w19, normalization.method = "LogNormalize", scale.factor = 10000)
w19 <- FindVariableGenes(object = w19, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = Inf, y.cutoff = 0.5)
length(x = w19@var.genes) #822
w19 <- ScaleData(object = w19, vars.to.regress = c("nUMI", "percentmito"))

#PCA/tSNE
w19 <- RunPCA(object = w19, pc.genes = w19@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = w19, pcs.use = 1:5)
PCHeatmap(object = w19, pc.use = 1, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = w19, pc.use = 2, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = w19, pc.use = 3, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
w19 <- JackStraw(object = w19, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = w19, PCs = 1:12) #1:8 significant
PCElbowPlot(object = w19)
w19 <- FindClusters(object = w19, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
w19 <- RunTSNE(object = w19, dims.use = 1:8, do.fast = T)
clustercolors = c("goldenrod","darkgreen","magenta")
TSNEPlot(object = w19, colors.use = clustercolors, pt.size = 1)
PCAPlot(object = w19, 1,2, cols.use = clustercolors)
PCAPlot(object = w19, 1,3, cols.use = clustercolors)
PCAPlot(object = w19, 3,2, cols.use = clustercolors)
PCAPlot(object = w19, 1,4, cols.use = clustercolors)
PCAPlot(object = w19, 1,5, cols.use = clustercolors)
PCAPlot(object = w19, 1,6, cols.use = clustercolors)
PCAPlot(object = w19, 1,7, cols.use = clustercolors)
PCAPlot(object = w19, 1,8, cols.use = clustercolors)
save(w19, file = "~/20180321-seuratw19.Robj")

#cell count and rename clusters
advgc_ids = WhichCells(w19, ident = "Adv. GC") #advanced gc, 88
length(w19@ident[(w19@ident %in% c('1'))]) #early gc, 60
length(w19@ident[(w19@ident %in% c('2'))]) #somatic cells, 39
w19@ident <- plyr::mapvalues(x = w19@ident, from = levels(w19@ident), to = c("Adv. GC", "Early GC", "Soma"))
TSNEPlot(object = w19, colors.use = clustercolors, pt.size = 1)

#gene expression datat
library(viridis)
scattercolors = c("cyan", "navy")
makegraphs = function(genelist){
  FeaturePlot(object = w19, features.plot = c(genelist), reduction.use = "pca", cols.use = scattercolors, no.legend = FALSE)
  FeaturePlot(object = w19, features.plot = c(genelist), reduction.use = "tsne", cols.use = scattercolors, no.legend = FALSE)
  VlnPlot(object = w19, features.plot = c(genelist), cols.use = clustercolors, x.lab.rot = T, size.x.use = 0)
}
makegraphs(c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT")) #germcellgenes
makegraphs(c("DNMT3A","DNMT3B","DNMT3L","MECP2","PRMT5")) #dnamethylgenes
makegraphs(c("EZH2", "EZH1", "SUZ12", "EED", "RBBP7", "RBBP4")) #Enzymes involved in H3K27meth
makegraphs(c("EHMT2","EHMT1","SETDB1","SUV39H1","SUV39H2")) #H3K9meth
makegraphs(c("NANOS2","LEFTY1","LEFTY2","PITX2","TEX14")) #male diff genes / dan's genes. RHOX6 and RHOX9 not found
makegraphs(c("PIWIL2", "HSP90AA1", "PIWIL4","MAEL","DDX4","MKI67")) #pirnacomps

save(w19, file = "20180322-seuratw19.Robj")

# read in the repeatome week 19 files ####
#figure out which cells are eliminated for poor quality
w19mt <- CreateSeuratObject(raw.data = txn19W_matrix, project = "mtpoor")
mitogenes <- grep(pattern = "^MT-", x = rownames(x = w19mt@data), value = TRUE)
percentmito <- Matrix::colSums(w19mt@raw.data[mitogenes, ])/Matrix::colSums(w19mt@raw.data)
w19mt <- AddMetaData(object = w19mt, metadata = percentmito, col.name = "percentmito")
length(x = w19@cell.names) #196
levels(w19mt@ident) = c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0")
w19mt = FilterCells(object = w19mt, subset.names = c("nGene", "percentmito"), low.thresholds = c(-Inf, .060), high.thresholds = c(Inf, 0.09))
GenePlot(object = w19mt, gene1 = "nUMI", gene2 = "percentmito", col.use = c("red", "magenta", "darkgreen", "cyan"))
mitohighcells = as.character(as.list(w19mt@data@Dimnames[[2]])) # these are the ones to remove from the te_19W matrix
rm(w19mt)
mitohighcells

#read in te_19W files and append as an assay
te_19W = read.delim(file = "te_19W.txt", row.names = 1)
te_19W = te_19W[,-grep("CCTCCTGA_25|ACATTGGC_27|AGAGTCAA_29|CTAGCATA_29|ATCCTGTA_30|CAAGACTA_30|CACTTCGA_30|CCTCCTGA_30|GAATCTGA_30",colnames(te_19W))]
w19 <- SetAssayData(w19, assay.type = "TE", slot = "raw.data", new.data = te_19W)
w19 <- NormalizeData(w19, assay.type = "TE", normalization.method = "LogNormalize", scale.factor = 10000)
w19 <- ScaleData(w19, assay.type = "TE", display.progress = T)
makegraphs(c("L1HS", "PIWIL4", "L1PREC2", "SVA_A", "HERVK9-int", "FordPrefect")) #re_1
par(mfrow = c(2, 2))
GenePlot(w19, "L1HS", "PIWIL4", col.use = clustercolors) #all three clusters
GenePlot(w19, "L1HS", "PIWIL2", col.use = clustercolors) #all three clusters
GenePlot(w19, "PIWIL4", "PIWIL2", col.use = clustercolors) #all three clusters
GenePlot(w19, "L1HS", "DNMT3A", col.use = clustercolors) #all three clusters

makegraphs(c("L1MD", "ALR/Alpha", "L1P", "L1PA7", "L1PA3", "LTR12D")) #re_2
makegraphs(c("L1PA10", "L1PA11", "L1PA12", "L1PA13", "L1PA14", "L1PA15", "L1PA15-16", "L1PA16")) #youngest set 1
makegraphs(c("L1PA17", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A")) #youngest set 2
makegraphs(c("L1ME1", "L1ME2", "L1ME2z", "L1ME3", "L1ME3A", "L1ME3B", "L1ME3C", "L1ME3Cz", "L1ME3D", "L1ME3E", "L1ME3F", "L1ME3G", "L1ME4a", "L1ME4b")) #old pt1
makegraphs(c("L1MEa", "L1MEb", "L1MEc", "L1MEd", "L1MEf", 'L1MEg', 'L1MEg1', 'L1MEg2', 'L1MEh', 'L1MEi', 'L1MEj', "L1ME4c", "L1ME5")) # old pt2
makegraphs(c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F")) #sva family, high in SC analysis
ave(w19, file = "20180322-seuratw19_txn_te.Robj")

load(file = "../20180323-seuratw19_txn_te.Robj")

# after 3/27/18 discussion ####
#L1HS colocalization with genes of interest
GenePlot(w19, "L1HS", "HSP90AA1", col.use = clustercolors) #all three clusters
GenePlot(w19, "L1HS", "MKI67", col.use = clustercolors) #all three clusters
GenePlot(w19, "L1HS", "POU5F1", col.use = clustercolors) #all three clusters

#PIWIL4 colocalization with other transposons

advgc_ids = names(w19@ident[(w19@ident %in% c('Adv. GC'))])
earlygc_ids = names(w19@ident[(w19@ident %in% c('Early GC'))])
allgc_ids = names(w19@ident[(w19@ident %in% c('Early GC','Adv. GC'))])

advgc = SubsetData(w19, cells.use = advgc_ids)
earlygc = SubsetData(w19, cells.use = earlygc_ids)
allgc = SubsetData(w19, cells.use = allgc_ids)

gcgeneplot = function(gene1, gene2){
  par(mfrow = c(2, 2))
  GenePlot(advgc, gene1, gene2, col.use = c('goldenrod'))
  GenePlot(earlygc, gene1, gene2, col.use = c('darkgreen'))
  GenePlot(allgc, gene1, gene2, col.use = clustercolors)
  GenePlot(w19, gene1, gene2, col.use = clustercolors)
}

par(mfrow = c(2, 2))
gcgeneplot("L1HS","HSP90AA1")
gcgeneplot("L1HS","MKI67")
gcgeneplot("L1HS","POU5F1")

gcgeneplot("PIWIL4","L1HS")
gcgeneplot("PIWIL4","SVA_A")
gcgeneplot("PIWIL4","L1PA3")
gcgeneplot("PIWIL4","L1PA7")
gcgeneplot("PIWIL4","L1PA11")
gcgeneplot("PIWIL4","L1PA14") 

