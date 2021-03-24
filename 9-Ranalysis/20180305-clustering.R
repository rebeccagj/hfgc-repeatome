setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/UMIcountwide/") #set wd to directory with UMI_Tools tsv file outputs

library(Seurat)
library(dplyr)
library(Matrix)

#read in data
all19W = read.csv("all19Wcounts.csv", header = T, row.names = 1, stringsAsFactors = F) # reads in the .csv with headers and with gene names as row names
colnames(all19W) #double checking that 26:30 are included in this cat'd file; there are 196 cells
length(colnames(all19W)) #double checking that 26:30 are included in this cat'd file; there are 196 cells
W19 <- CreateSeuratObject(raw.data = all19W, min.cells = 3, min.genes = 200, project = "W19hFGC") #create Seurat object

#filter and normalize raw data
W19 <- FilterCells(W19, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
W19 <- NormalizeData(W19)
W19 <- ScaleData(W19, display.progress = T)
#mito.genes <- read.delim("/Users/Rebecca/Downloads/mito_genes.txt",header = F)
#mito.genes = as.vector(mito.genes$V1)
#percent.mito <- Matrix::colSums(W19@raw.data[mito.genes, ])/Matrix::colSums(W19@raw.data)
W19 <- AddMetaData(object = W19, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = W19, features.plot = c("nGene", "nUMI"), nCol = 2)

W19 <- NormalizeData(object = W19, normalization.method = "LogNormalize", scale.factor = 10000)

W19 <- FindVariableGenes(object = W19, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

W19 <- ScaleData(object = W19, vars.to.regress = c("nUMI"))
W19 <- RunPCA(object = W19, pc.genes = W19@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = W19, dim.1 = 1, dim.2 = 2)
PCHeatmap(object = W19, pc.use = 1, cells.use = 197, do.balanced = TRUE, label.columns = FALSE)
PCElbowPlot(object = W19)
W19 <- FindClusters(object = W19, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
W19 <- RunTSNE(object = W19, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = W19)
