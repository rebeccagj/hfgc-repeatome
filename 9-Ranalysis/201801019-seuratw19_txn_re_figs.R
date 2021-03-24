# better gtf analysis! relevant bits for figure are lines 137-260
setwd("~/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_RE-wk19/9-Ranalysis/") #set wd to directory with UMI_Tools tsv file outputs

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
read.delim(file = "txn_19W.txt")
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
read.delim("te_19W.txt")
#also i changed 1:1 from "rn" to "genes" in the txt file manually via a text editor

#week 19 transcriptome dataset input and clustering ####
setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/genebycell_19W/") #set wd to directory with UMI_Tools tsv file outputs
library(Seurat)
txn_te_W19_table = read.delim(file = "txn_te_19W.txt", row.names = 1) #re/txn combo files
txn_te_W19_matrix = as(as.matrix(txn_te_W19_table), "dgCMatrix") #makes a sparse matrix of txn+re dataframe

#qc filtering
txn_re_w19 <- CreateSeuratObject(raw.data = txn_te_W19_matrix, project = "hFGC19retxn")
mitogenes <- grep(pattern = "^MT-", x = rownames(x = txn_re_w19@data), value = TRUE)
percentmito <- Matrix::colSums(txn_re_w19@raw.data[mitogenes, ])/Matrix::colSums(txn_re_w19@raw.data)
txn_re_w19 <- AddMetaData(object = txn_re_w19, metadata = percentmito, col.name = "percentmito")
length(x = txn_re_w19@cell.names) #196
levels(txn_re_w19@ident) = c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0")
par(mfrow = c(1, 2))
GenePlot(object = txn_re_w19, gene1 = "nUMI", gene2 = "percentmito", col.use = c("red", "magenta", "darkgreen", "cyan"))
GenePlot(object = txn_re_w19, gene1 = "nUMI", gene2 = "nGene")
txn_re_w19 = FilterCells(object = txn_re_w19, subset.names = c("nGene", "percentmito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(Inf, 0.06))
txn_re_w19 #An object of class seurat in project hFGC19. 25657 genes across 191 samples.

#normalization and variable genes
txn_re_w19 <- NormalizeData(object = txn_re_w19, normalization.method = "LogNormalize", scale.factor = 10000)
par(mfrow = c(1, 1))
txn_re_w19 <- FindVariableGenes(object = txn_re_w19, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = Inf, y.cutoff = 0.5)
length(x = txn_re_w19@var.genes) #782
txn_re_w19 <- ScaleData(object = txn_re_w19, vars.to.regress = c("nUMI", "percentmito"))

#PCA/tSNE
txn_re_w19 <- RunPCA(object = txn_re_w19, pc.genes = txn_re_w19@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 30)
VizPCA(object = txn_re_w19, pcs.use = 1:5)
PCHeatmap(object = txn_re_w19, pc.use = 1, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = txn_re_w19, pc.use = 2, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = txn_re_w19, pc.use = 3, cells.use = 196, do.balanced = TRUE, label.columns = FALSE)
txn_re_w19 <- JackStraw(object = txn_re_w19, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = txn_re_w19, PCs = 1:12) #1:10 significant
PCElbowPlot(object = txn_re_w19)
txn_re_w19 <- FindClusters(object = txn_re_w19, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
txn_re_w19 <- RunTSNE(object = txn_re_w19, dims.use = 1:10, do.fast = T)
clustercolors = c("goldenrod","darkgreen","magenta")
TSNEPlot(object = txn_re_w19, colors.use = clustercolors, pt.size = 1)
PCAPlot(object = txn_re_w19, 1,2, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,3, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 3,2, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,4, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,5, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,6, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,7, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,8, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,9, cols.use = clustercolors)
PCAPlot(object = txn_re_w19, 1,10, cols.use = clustercolors)
length(txn_re_w19@ident[(txn_re_w19@ident %in% c('0'))]) #adv gc, 89
length(txn_re_w19@ident[(txn_re_w19@ident %in% c('1'))]) #early gc, 63
length(txn_re_w19@ident[(txn_re_w19@ident %in% c('2'))]) #somatic cells, 39
txn_re_w19@ident <- plyr::mapvalues(x = txn_re_w19@ident, from = levels(txn_re_w19@ident), to = c("Adv. GC", "Early GC", "Soma"))
txn_re_w19 <- RunTSNE(object = txn_re_w19, dims.use = 1:10, do.fast = T)
TSNEPlot(object = txn_re_w19, colors.use = clustercolors, pt.size = 1)
save(txn_re_w19, file = "../20180420-seuratw19_txn_re_figs.Robj")

############################################################################
#most recent iteration for generating figures, contains both re and txn data
############################################################################
setwd("~/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_RE-wk19/9-Ranalysis/")
load(file = "20180420-seuratw19_txn_re_figs.Robj") 

scattercolors = c("cyan", "navy")
clustercolors = c("blue","forestgreen","gray45")
TSNEPlot(object = txn_re_w19, colors.use = clustercolors, pt.size = 2)

advgc_so = SubsetData(txn_re_w19, cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))]))
earlygc_so = SubsetData(txn_re_w19, cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))]))
allgc_so = SubsetData(txn_re_w19, cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC','Adv. GC'))]))

#gene expression data
makegraphs = function(genelist){
  #FeaturePlot(object = txn_re_w19, features.plot = c(genelist), reduction.use = "pca", cols.use = scattercolors, no.legend = FALSE, pt.size=2)
  par(mfrow = c(2, 3))
  FeaturePlot(object = txn_re_w19, features.plot = c(genelist), reduction.use = "tsne", cols.use = scattercolors, no.legend = FALSE, pt.size=2)
  #VlnPlot(object = txn_re_w19, features.plot = c(genelist), cols.use = clustercolors, x.lab.rot = T, size.x.use = 0, point.size.use = 2, return.plotlist = TRUE, same.y.lims = T)
}
gcgeneplot = function(gene1, gene2){
  par(mfrow = c(2, 2))
  GenePlot(advgc_so, gene1, gene2, col.use = c('blue'), cex.use = 1.5)
  GenePlot(earlygc_so, gene1, gene2, col.use = c('forestgreen'),cex.use = 1.5)
  GenePlot(allgc_so, gene1, gene2, col.use = clustercolors, cex.use = 1.5)
  GenePlot(txn_re_w19, gene1, gene2, col.use = clustercolors, cex.use = 1.5)
}
makegraphs(c("DNMT3A","DNMT3B","DNMT3L","MECP2","PRMT5")) #dnamethylgenes
makegraphs(c("EZH2", "EZH1", "SUZ12", "EED", "RBBP7", "RBBP4")) #Enzymes involved in H3K27meth
makegraphs(c("EHMT2","EHMT1","SETDB1","SUV39H1","SUV39H2")) #H3K9meth
makegraphs(c("NANOS2","LEFTY1","LEFTY2","PITX2","TEX14")) #male diff genes / dan's genes. RHOX6 and RHOX9 not found
makegraphs(c("L1HS", "PIWIL4", "L1PREC2", "SVA_A", "HERVK9-int", "FordPrefect")) 
makegraphs(c("L1MD", "ALR/Alpha", "L1P", "L1PA7", "L1PA3", "LTR12D")) #re_2
makegraphs(c("L1PA10", "L1PA11", "L1PA12", "L1PA13", "L1PA14", "L1PA15", "L1PA15-16", "L1PA16")) #youngest set 1
makegraphs(c("L1PA17", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA8A")) #youngest set 2
makegraphs(c("L1ME1", "L1ME2", "L1ME2z", "L1ME3", "L1ME3A", "L1ME3B", "L1ME3C", "L1ME3Cz", "L1ME3D", "L1ME3E", "L1ME3F", "L1ME3G", "L1ME4a", "L1ME4b")) #old pt1
makegraphs(c("L1MEa", "L1MEb", "L1MEc", "L1MEd", "L1MEf", 'L1MEg', 'L1MEg1', 'L1MEg2', 'L1MEh', 'L1MEi', 'L1MEj', "L1ME4c", "L1ME5")) # old pt2
makegraphs(c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F")) #sva family, high in SC analysis
makegraphs(c("MAGEB2","RHOXF2B","KIF3A","RHOXF2","ID1","PABPC4")) #bioinformatically determined markers of adv. GC cluster

#for pdfs, pt size 2
makegraphs(c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT","AMH","WT1","SOX9","ARX","TCF21","CYP17A1")) #germcellgenes and somatic genes for same size export
makegraphs(c("POU5F1", "DDX4", "MAEL", "HSP90AA1", "PIWIL2", "PIWIL4")) #germcellgenes
makegraphs(c("AMH","WT1","SOX9","ARX","TCF21","CYP17A1")) #sertoli markers from LiLi etal
makegraphs(c("ARX","TCF21","CYP17A1")) #leydig markers from LiLi etal
makegraphs(c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT")) #germcellgenes
makegraphs(c("L1PA2", "L1HS", "L1PA3","DAZL")) #re genes
makegraphs(c("PIWIL2", "PIWIL4","MAEL","HSP90AA1")) #pirnacomps
gcgeneplot("L1HS","HSP90AA1")
gcgeneplot("L1HS","PIWIL2")
gcgeneplot("L1HS","PIWIL4")
gcgeneplot("HSP90AA1","PIWIL4")

#for sig values of the gene expression plots ####
sigvaluesfxn = function(gene1, gene2){
  advgcdata = FetchData(advgc_so, vars.all = c(gene1,gene2))
  write.table(advgcdata, file = "advgcdata.txt", quote = FALSE, sep = "\t", row.names = TRUE)
  
  earlygcdata = FetchData(earlygc_so, vars.all = c(gene1,gene2))
  write.table(earlygcdata, file = "earlygcdata.txt", quote = FALSE, sep = "\t", row.names = TRUE)
  
  allgcdata = FetchData(allgc_so, vars.all = c(gene1,gene2))
  write.table(allgcdata, file = "allgcdata.txt", quote = FALSE, sep = "\t", row.names = TRUE)
  
  somaandgc = FetchData(txn_re_w19, vars.all = c(gene1,gene2))
  write.table(somaandgc, file = "somaandgc.txt", quote = FALSE, sep = "\t", row.names = TRUE)
}
setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/genebycell_19W/")
setwd("gene_coor_stats/l1hs_hsp90a/")
sigvaluesfxn("L1HS","HSP90AA1")
setwd("../l1hs_piwil2/")
sigvaluesfxn("L1HS","PIWIL2")
setwd("../l1hs_piwil4/")
sigvaluesfxn("L1HS","PIWIL4")

# read in the txn alone and repeatome alone week 19 files ####
#figure out which cells are eliminated for poor quality
txn_w19_mtcheck <- CreateSeuratObject(raw.data = txn_W19_matrix, project = "mtpoor")
mitogenes <- grep(pattern = "^MT-", x = rownames(x = txn_w19_mtcheck@data), value = TRUE)
percentmito <- Matrix::colSums(txn_w19_mtcheck@raw.data[mitogenes, ])/Matrix::colSums(txn_w19_mtcheck@raw.data)
txn_w19_mtcheck <- AddMetaData(object = txn_w19_mtcheck, metadata = percentmito, col.name = "percentmito")
length(x = txn_w19_mtcheck@cell.names) #196
levels(txn_w19_mtcheck@ident) = c("0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0")
txn_w19_mtcheck = FilterCells(object = txn_w19_mtcheck, subset.names = c("nGene", "percentmito"), low.thresholds = c(-Inf, .060), high.thresholds = c(Inf, 0.09))
GenePlot(object = txn_w19_mtcheck, gene1 = "nUMI", gene2 = "percentmito", col.use = c("red", "magenta", "darkgreen", "cyan"))
mitohighcells = as.character(as.list(txn_w19_mtcheck@data@Dimnames[[2]])) # these are the ones to remove from the te_19W matrix
rm(txn_w19_mtcheck)
mitohighcells #next line greps out the 'mitohighcells' string

#read in re_W19 alone files
re_W19_table = read.delim(file = "te_19W.txt", row.names = 1) #re alone files
re_W19_matrix = as(as.matrix(re_W19_table), "dgCMatrix") #makes a sparse matrix of re only dataframe
re_W19_matrix = re_W19_matrix[,-grep("CCTCCTGA_25|ACATTGGC_27|AGAGTCAA_29|CTAGCATA_29|ATCCTGTA_30|CAAGACTA_30|CACTTCGA_30|CCTCCTGA_30|GAATCTGA_30",colnames(re_W19_matrix))] #greps out the 'mitohighcells' string from re matrix

#read in the txn_W19 alone files
txn_W19_table = read.delim(file = "txn_19W.txt", row.names = 1) #txn alone files
txn_W19_matrix = as(as.matrix(txn_W19_table), "dgCMatrix") #makes a sparse matrix of txn only dataframe
txn_W19_matrix = txn_W19_matrix[,-grep("CCTCCTGA_25|ACATTGGC_27|AGAGTCAA_29|CTAGCATA_29|ATCCTGTA_30|CAAGACTA_30|CACTTCGA_30|CCTCCTGA_30|GAATCTGA_30",colnames(txn_W19_matrix))]  #greps out the 'mitohighcells' string from txn matrix

#append the re_W19 file as an assay of a txn_W19 seurat object
txn_assayRE_w19 <- CreateSeuratObject(raw.data = txn_W19_matrix, project = "hFGC19txnassay")
txn_assayRE_w19 <- NormalizeData(object = txn_assayRE_w19, normalization.method = "LogNormalize", scale.factor = 10000)
txn_assayRE_w19 <- FindVariableGenes(object = txn_assayRE_w19, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = Inf, y.cutoff = 0.5)
length(x = txn_assayRE_w19@var.genes) #782
txn_assayRE_w19 <- ScaleData(object = txn_assayRE_w19, vars.to.regress = c("nUMI"))
txn_assayRE_w19 <- SetAssayData(txn_assayRE_w19, assay.type = "RE", slot = "raw.data", new.data = re_W19_matrix)
txn_assayRE_w19 <- NormalizeData(txn_assayRE_w19, assay.type = "RE", normalization.method = "LogNormalize", scale.factor = 10000)
txn_assayRE_w19 <- ScaleData(txn_assayRE_w19, assay.type = "RE", display.progress = T)
txn_assayRE_w19 <- RunPCA(object = txn_assayRE_w19, pc.genes = txn_assayRE_w19@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 30)
txn_assayRE_w19 <- JackStraw(object = txn_assayRE_w19, num.replicate = 100, do.print = FALSE)
txn_assayRE_w19 <- FindClusters(object = txn_assayRE_w19, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
txn_assayRE_w19 <- RunTSNE(object = txn_assayRE_w19, dims.use = 1:10, do.fast = T)
clustercolors = c("goldenrod","darkgreen","magenta")
TSNEPlot(object = txn_assayRE_w19, colors.use = clustercolors, pt.size = 1)

makegraphs_assay = function(genelist){
  FeaturePlot(object = txn_assayRE_w19, features.plot = c(genelist), reduction.use = "pca", cols.use = scattercolors, no.legend = FALSE)
  FeaturePlot(object = txn_assayRE_w19, features.plot = c(genelist), reduction.use = "tsne", cols.use = scattercolors, no.legend = FALSE)
  VlnPlot(object = txn_assayRE_w19, features.plot = c(genelist), cols.use = clustercolors, x.lab.rot = T, size.x.use = 0)
}
makegraphs_assay(c("L1HS", "PIWIL4", "L1PREC2", "SVA_A", "HERVK9-int", "FordPrefect")) #seems similar to the original object with everything pre-included, going back to that one


#PIWIL4 and L1HS colocalization with other transposons ####
#advgc_ids_re = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))])
#earlygc_ids_re = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))])
#allgc_ids_re = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC','Adv. GC'))])

gcgeneplot("L1HS","HSP90AA1")
gcgeneplot("L1HS","PIWIL2")
gcgeneplot("L1HS","PIWIL4")
gcgeneplot("HSP90AA1","PIWIL4")
gcgeneplot("PIWIL4","HSP90AA1")

gcgeneplot("L1PA2","PIWIL2")
gcgeneplot("L1PA2","PIWIL4")
gcgeneplot("L1PA3","PIWIL2")
gcgeneplot("L1PA3","PIWIL4")

gcgeneplot("L1PA2","HSP90AA1") #checking to see if L1HS or these are better, decided not to follow this up bs L1HS neg coor is better
gcgeneplot("L1PA3","HSP90AA1")

#pre-figure investigations ####
gcgeneplot("L1HS","MKI67")
gcgeneplot("L1HS","POU5F1")
gcgeneplot("PIWIL4","L1HS")
gcgeneplot("PIWIL4","SVA_A")
gcgeneplot("PIWIL4","L1PA3")
gcgeneplot("PIWIL4","L1PA7")
gcgeneplot("PIWIL4","L1PA11")
gcgeneplot("PIWIL4","L1PA14") 

gcgeneplot("PIWIL4","AluY")
gcgeneplot("PIWIL4","AluYa5")
gcgeneplot("PIWIL4","AluYa8")
gcgeneplot("PIWIL4","AluYb8")
gcgeneplot("PIWIL4","AluYb9")
gcgeneplot("PIWIL4","AluYc")
gcgeneplot("PIWIL4","AluYc3")
gcgeneplot("PIWIL4","AluYd8")
gcgeneplot("PIWIL4","AluYe5")
gcgeneplot("PIWIL4","AluYe6")
gcgeneplot("PIWIL4","AluYf1")
gcgeneplot("PIWIL4","AluYg6")
gcgeneplot("PIWIL4","AluYh3")
gcgeneplot("PIWIL4","AluYh3a3")
gcgeneplot("PIWIL4","AluYh7")
gcgeneplot("PIWIL4","AluYh9")
gcgeneplot("PIWIL4","AluYi6")
gcgeneplot("PIWIL4","AluYi6_4d")
gcgeneplot("PIWIL4","AluYj4")
gcgeneplot("PIWIL4","AluYk11")
gcgeneplot("PIWIL4","AluYk12")
gcgeneplot("PIWIL4","AluYk2")
gcgeneplot("PIWIL4","AluYk3")
gcgeneplot("PIWIL4","AluYk4")
gcgeneplot("PIWIL4","AluYm1")
gcgeneplot("L1HS","DDX4")
gcgeneplot("HIWI2","L1HS")

# 2018 12 04 resubmission request

sigvaluesfxn("L1PA2","L1PA3")
setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_RE-wk19/9-Ranalysis/genebycell_19W/gene_coor_stats/l1pa2_l1pa3/")
sigvaluesfxn("L1PA2","L1PA3")

# all germ cell analysis - nothing different####
allgc_re <- ScaleData(object = allgc_re, vars.to.regress = c("nUMI"))
allgc_re <- RunPCA(object = allgc_re, pc.genes = allgc_re@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PCHeatmap(object = allgc_re, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)
allgc_re <- JackStraw(object = allgc_re, num.replicate = 100)
JackStrawPlot(object = allgc_re, PCs = 1:12)
PCElbowPlot(allgc_re)
allgc_re <- FindClusters(object = allgc_re, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
TSNEPlot(object = allgc_re)

#ugh nothing is different. going back to original R object txn_re_w19 so that all the layouts match :'( 
# differential expression between txn and re in txn_re_w19####  
#20180420 update: Future Rebecca please append RE and process it to txn_re_w19
#or can just reload the 3/28files

TSNEPlot(object = txn_re_w19, colors.use = clustercolors, pt.size = 1)
advGCmarkers <- FindMarkers(object = txn_re_w19, ident.1 = "Adv. GC", min.pct = 0.25, test.use = "roc")
earlyGCmarkers <- FindMarkers(object = txn_re_w19, ident.1 = "Early GC", min.pct = 0.25, test.use = "roc")
re_advGCmarkers <- FindMarkers(object = txn_re_w19, assay.type = "RE", ident.1 = "Adv. GC", min.pct = 0.25, test.use = "roc")
re_earlyGCmarkers <- FindMarkers(object = txn_re_w19, assay.type = "RE", ident.1 = "Early GC", min.pct = 0.25, test.use = "roc")
write.table(advGCmarkers, file = "20180328-advGCmarkers.txt", quote = FALSE, sep = "\t")
write.table(earlyGCmarkers, file = "20180328-earlyGCmarkers.txt", quote = FALSE, sep = "\t")
write.table(re_advGCmarkers, file = "20180328-re_advGCmarkers.txt", quote = FALSE, sep = "\t")
write.table(re_earlyGCmarkers, file = "20180328-re_earlyGCmarkers.txt", quote = FALSE, sep = "\t")

#read these files in if starting from a fresh environment
advGCmarkers = read.delim(file = "20180328-advGCmarkers.txt", sep = "\t", row.names = 1) #adv gc markers from roc test
earlyGCmarkers = read.delim(file = "20180328-earlyGCmarkers.txt", sep = "\t", row.names = 1) #adv gc markers from roc test
re_advGCmarkers = read.delim("20180328-re_advGCmarkers.txt", row.names = 1)
re_earlyGCmarkers = read.delim("20180328-re_earlyGCmarkers.txt", row.names = 1)

#top 25 differentially expressed RE from the adv GC population ####
re_advGCmarkers <- re_advGCmarkers[- grep("tRNA*", rownames(re_advGCmarkers)),] #removed tRNAs to get all non-tRNA genes
re_advGCmarkers = tibble::rownames_to_column(re_advGCmarkers, "Gene") #appends genes from rowname to column 1, otherwise next line will discard them
library(dplyr)
re_advGC_top25roc_hm = re_advGCmarkers %>% top_n(25, re_advGCmarkers$myAUC) #grabs the top 25 DE REs by myAUC sort
re_advGC_top25roc_hm = re_advGC_top25roc_hm$Gene

txn_hmp_genes = c("POU5F1",'DDX4',"DAZL","NANOS2","PIWIL1","PIWIL2",'PIWIL3','PIWIL4','MAEL','HSP90AA1','TDRD1','TDRD5','TDRD6','TDRD9','TDRD12','TDRKH','MORC1','PLD6','HENMT1','RNF17')
heatmap_figure_genes = txn_hmp_genes
heatmap_figure_genes[22:46] = re_advGC_top25roc_hm
heatmap_figure_genes

DoHeatmap(txn_re_w19, genes.use = heatmap_figure_genes, slim.col.label = TRUE,
          remove.key = FALSE, group.label.rot = F, rotate.key = F, cex.row = 8,
          col.low = "black", col.mid = "blue", col.high = "yellow", group.spacing = 0.5)

DoHeatmap(txn_re_w19, genes.use = re_advGC_top25roc_hm, slim.col.label = TRUE,
          remove.key = FALSE, group.label.rot = F, rotate.key = F, cex.row = 8,
          col.low = "black", col.mid = "blue", col.high = "yellow", group.spacing = 0.5)

DoHeatmap(txn_re_w19, genes.use = txn_hmp_genes, slim.col.label = TRUE,
          remove.key = FALSE, group.label.rot = F, rotate.key = F, cex.row = 8,
          col.low = "black", col.mid = "blue", col.high = "yellow", group.spacing = 0.5)

# Heatmap of L1PAs, gene lists from  doi: 10.1186/s13100-017-0107-y ####

L1Pgenes = c("L1HS", "L1P4a_5end", "L1P4b_5end", "L1P4c_5end", "L1P4d_5end", "L1P4e_5end","L1PA2", 'L1PA3',
             "L1PA4", "L1PA5", "L1PA6", "L1PA7", "L1PA8", "L1PA10", "L1PA11", "L1PA12", "L1PA12_5", "L1PA13", "L1PA13_5",
             "L1PA14", "L1PA14_5", "L1PA15", "L1PA16", "L1PA16_5", "L1PA17_5", "L1PA7_5",  "L1 PB1", "L1 PB2", "L1PB2c", 
             "L1 PB3", "L1 PB4", "L1PBA1_5", "L1PBA_5", "L1PBB_5", "L1PREC1", "L1PREC2", "L1P_MA2", "L1")

DoHeatmap(txn_re_w19, genes.use = L1Pgenes, slim.col.label = TRUE,
          remove.key = FALSE, group.label.rot = F, rotate.key = F, cex.row = 8,
          col.low = "black", col.mid = "blue", col.high = "yellow", group.spacing = 0.5)

# miRNA processing/effector pathway genes from Boris in response to reviewers ####

class(miRNApathwaygenes) = c("DICER1", "DROSHA", "DGCR8", "XPO5", "R3HDM1", "AGO1", "AGO2", "AGO3", "AGO4",  "ARPP21")

VlnPlot(txn_re_w19, miRNApathwaygenes, cols.use = clustercolors)

DoHeatmap(txn_re_w19, genes.use = miRNApathwaygenes, slim.col.label = TRUE,
          remove.key = FALSE, group.label.rot = F, rotate.key = F, cex.row = 8,
          col.low = "black", col.mid = "blue", col.high = "yellow", group.spacing = 0.5)

gcgeneplot("L1HS",miRNApathwaygenes[1])
gcgeneplot("L1HS",miRNApathwaygenes[2])
gcgeneplot("L1HS",miRNApathwaygenes[3])
gcgeneplot("L1HS",miRNApathwaygenes[4])
gcgeneplot("L1HS",miRNApathwaygenes[5])
gcgeneplot("L1HS",miRNApathwaygenes[6])
gcgeneplot("L1HS",miRNApathwaygenes[7])
gcgeneplot("L1HS",miRNApathwaygenes[8])
gcgeneplot("L1HS",miRNApathwaygenes[9])
gcgeneplot("L1HS",miRNApathwaygenes[10])

gcgeneplot("PIWIL4",miRNApathwaygenes[1])
gcgeneplot("PIWIL4",miRNApathwaygenes[2])
gcgeneplot("PIWIL4",miRNApathwaygenes[3])
gcgeneplot("PIWIL4",miRNApathwaygenes[4])
gcgeneplot("PIWIL4",miRNApathwaygenes[5])
gcgeneplot("PIWIL4",miRNApathwaygenes[6])
gcgeneplot("PIWIL4",miRNApathwaygenes[7])
gcgeneplot("PIWIL4",miRNApathwaygenes[8])
gcgeneplot("PIWIL4",miRNApathwaygenes[9])
gcgeneplot("PIWIL4",miRNApathwaygenes[10])

gcgeneplot("PIWIL4","DDX4")
z
ciliagenes = read.table("~/Downloads/cilia-hmn-genes.txt")
cilialist = as.character(ciliagenes$V1)

ciliashort = c('ARL13B', 'BBS4', 'IFT88', 'SEPT7', 'RAB8A')

VlnPlot(txn_re_w19, ciliashort, cols.use = clustercolors)
FeaturePlot(txn_re_w19, ciliashort, cols.use = scattercolors, no.legend = F)

ciliahomo = read.table("~/Downloads/mapped_orthologs.txt", fill = T)
cilialist = as.character(ciliahomo$V22)

# same size dots and WxH for new figure
quaantiles = draw_quantiles = c(0.25, 0.5, 0.75)
VlnPlot(txn_re_w19, c("L1PA2"), cols.use = clustercolors, point.size.use = 0.5)
VlnPlot(txn_re_w19, c("L1HS"), cols.use = clustercolors, y.max = 1.6, point.size.use = 0.5, do.return = T)
VlnPlot(txn_re_w19, c("L1PA3"), cols.use = clustercolors, point.size.use = 0.5)

install.packages('gridExtra')
library(gridExtra)
p1 = DBPlot("L1HS", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,1.6,0.4), limits = c(0,1.6))
p2 = DBPlot("L1PA2", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,2.5,0.5), limits = c(0,2.5))
p3 = DBPlot("L1PA3", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,2.5,0.5), limits = c(0,2.5))
p4 = DBPlot("HSP90AA1", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Adv. GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,4.25,1), limits = c(0,4.25))
grid.arrange(p1, p2, p3, p4, nrow = 1)
p5 = DBPlot("L1HS", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,2,0.5), limits = c(0,2.1))
p6 = DBPlot("L1PA2", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,2.3,0.5), limits = c(0,2.3))
p7 = DBPlot("L1PA3", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,2.5,0.5), limits = c(0,2.5))
p8 = DBPlot("HSP90AA1", "txn_re_w19", group.by = "ident", color.by = "ident", cells.use = names(txn_re_w19@ident[(txn_re_w19@ident %in% c('Early GC'))]),
       plots = c("vlnplot","jitter"), size = 1) + scale_y_continuous(breaks=seq(0,4.25,1), limits = c(0,4.25))
grid.arrange(p5, p6, p7, p8, nrow = 1)
