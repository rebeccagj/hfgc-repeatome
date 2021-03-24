setwd("/Users/Rebecca/Documents/Graduate_School/Laird_Lab/Experiments/RJ31-LiLi_et_al_2017-HFGC_TE-alignment/9-Ranalysis/GSE86146_RAW/") #set wd to directory with UMI_Tools tsv file outputs

#function for inputting data and assigning to global environment ####
readrenamecol = function(arg1, arg2){
  arg1 = read.delim(file = arg2, header = T, sep = "\t", row.names = 1)
} 

#read in the 6, 19W samples ####
gsm28 = readrenamecol(gsm28, "GSM2306028_M_19W_embryo1_101_gene_expression.txt")
gsm29 = readrenamecol(gsm29, "GSM2306029_M_19W_embryo1_24_gene_expression.txt")
gsm30 = readrenamecol(gsm30, "GSM2306030_M_19W_embryo1_26_gene_expression.txt")
gsm31 = readrenamecol(gsm31, "GSM2306031_M_19W_embryo2_102_gene_expression.txt")
gsm32 = readrenamecol(gsm32, "GSM2306032_M_19W_embryo2_103_gene_expression.txt")
gsm33 = readrenamecol(gsm33, "GSM2306033_M_19W_embryo2_104_gene_expression.txt")
class(gsm28) #data.frame

#concatinate all datasets, dropping rows that are not present in all SRR sets ####
all19W = merge.data.frame(gsm28, gsm29)
all19W = merge.data.frame(all19W, gsm30)
all19W = merge.data.frame(all19W, gsm31)
all19W = merge.data.frame(all19W, gsm32)
all19W = merge.data.frame(all19W, gsm33)
write.table(all19W, file = "GSE86146_RAW_all19Wcounts.txt", sep = "\t", quote = F, col.names = NA) #write merged data.frame to csv file
rm(gsm28,gsm29,gsm30,gsm31,gsm32,gsm33)

# start analysis ####
library(Seurat)

all19W = read.csv("GSE86146_RAW_all19Wcounts.txt", sep = "\t", row.names = 2)
all19W = all19W[-1]
all19W_matrix = as(as.matrix(all19W), "dgCMatrix")
w19 <- CreateSeuratObject(raw.data = all19W_matrix, min.cells = 5, is.expr = 1, project = "hFGC19")
VlnPlot(object = w19, features.plot = c("nGene", "nUMI"), nCol = 2)
GenePlot(object = w19, gene1 = "nUMI", gene2 = "nGene")
w19 <- NormalizeData(object = w19, normalization.method = "LogNormalize", scale.factor = 10000)
w19 <- FindVariableGenes(object = w19, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = Inf, y.cutoff = 0.5)
length(x = w19@var.genes) #587
w19 <- ScaleData(object = w19, vars.to.regress = c("nUMI"))
w19 #17344 genes across 195 samples
w19 <- RunPCA(object = w19, pc.genes = w19@var.genes, pcs.compute = 20, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = w19, dim.1 = 1, dim.2 = 3)
PCHeatmap(object = w19, pc.use = 1:5, do.balanced = TRUE, label.columns = F, remove.key = FALSE )
PCHeatmap(object = w19, pc.use = 1, do.balanced = TRUE, label.columns = F, remove.key = FALSE, col.use = PurpleAndYellow())
w19 <- JackStraw(object = w19, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = w19, PCs = 1:12) #seems 1-10 are sig
PCElbowPlot(object = w19)
w19 <- FindClusters(object = w19, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
w19 <- RunTSNE(object = w19, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(w19, do.label = T)
VlnPlot(object = w19, features.plot = c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT"))
FeaturePlot(object = w19, features.plot = c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT"), reduction.use = "tsne")
FeaturePlot(object = w19, features.plot = c("NANOG", "DPPA3", "POU5F1", "DAZL", "DDX4", "KIT"), reduction.use = "pca")
FeaturePlot(object = w19, features.plot = c("MKI67", "HSP90AA1", "PIWIL4"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = w19, features.plot = c("MKI67", "HSP90AA1", "PIWIL4"), reduction.use = "pca", no.legend = FALSE)
VlnPlot(object = w19, features.plot = c("MKI67", "HSP90AA1", "PIWIL4"))


cluster0 = names(w19@ident[(w19@ident %in% c('0'))])
cluster1 = names(w19@ident[(w19@ident %in% c('1'))])
cluster2 = names(w19@ident[(w19@ident %in% c('2'))])
