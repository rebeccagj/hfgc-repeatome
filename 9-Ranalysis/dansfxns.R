#These are my (sc)RNAseq functions and some things that they rely on.

MYcolors <- c("forestgreen", "#FFFFFF", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#FFFFFF", "#C38700", "#4999C6", "#008662")

library(ggplot2, quietly = T)
library(dplyr, quietly = T)

DBPlot <- function(var, object = DEFAULT, main = NULL, sub = NULL, group.by = "Tage",
                   color.by = "age", cells.use = NULL,
                   size=0.1, shape=16, low = "#F0E442", high = "#0072B2", colors = c(1:8),
                   ylab = NULL, xlab = NULL, plots = c("jitter","vlnplot"), hline=NULL, labels = NULL,
                   range = NULL){
  #Makes a plot where color is overlayed onto the dim.reduction plot of choice.
  #
  #object                 the Seurat Object
  #var                    Target Variable = either the metadata name or values 
  #group.by               metadata to use for separating values.  Default is by sample.
  #color.by               metadata to use for coloring.
  #cells.use              cells to include
  #main                   plot title
  #sub                    plot subtitle
  #colors                 indexes of color from MYcolors
  #plots                  types of plots to include: possibilities = "jitter", "boxplot", "vlnplot"
  #labels                 names to change x labels to.  Defaults to my current samples
  #ylab                   y axis label, default is 
  #hline                  value(s) where a dashed horizontal line should go
  
  #If cells.use = NA (was not provided), populate it to be all cells or all samples.
  if (classof(object)=="seurat" & is.null(cells.use)) {cells.use <- eval(expr = parse(text = paste0(object,"@cell.names")))}
  if (classof(object)=="RNAseq" & is.null(cells.use)) {cells.use <- meta("Samples")}
  
  ###Determine what the y-axis should be.
  #non-direct input options are the name of a metadata or a gene.
  #Both these options also get a default axis title.
  if (is.meta(var, object)){Y <- meta(var, object)}
  if (is.gene(var, object)){Y <- gene(var, object)}
  if (!(is.meta(var, object)|is.gene(var, object))){Y <- var}
  #Add default y-labels if the var is a metadata or a gene name.
  if (is.null(ylab) & is.meta(var, object)) {ylab <- var}
  if (is.null(ylab) & is.gene(var,object)) {ylab <- paste0(var," expression")}
  ###Make dataframe for storing the plotting data:
  #The data =Y, how to group the data = sample, how to color the groupings = color, and what the shape should be if there is a "jitter" made.
  full_dat <- data.frame(Y = Y,
                         sample = meta(group.by, object),
                         color = meta(color.by, object),
                         Shape = shape)
  #Subset the data.frame to only the cells in cell.use.
  if (classof(object)=="seurat"){
    Target_dat <- full_dat[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use,]
  }
  if (classof(object)=="RNAseq"){
    Target_dat <- full_dat[meta("Samples") %in% cells.use,]
  }
  #####Start making the plot
  p <- ggplot(Target_dat, aes(x=sample, y=Y, fill= color, color = color)) +
    #And set the legend to not have a title
    theme(legend.title=element_blank())
  #Add the y label
  p<- p + ylab(ylab)
  #Set the y-axis limits if a range is given.
  if (!is.null(range)){p <- p + ylim(range)}
  ###Add data based on what is requested in plots, ordered by their order
  for (i in 1:length(plots)){
    if (plots[i] == "boxplot") {p <- p + geom_boxplot(width=0.2, aes(color = "white"), outlier.shape = NA)} 
    if (plots[i] == "jitter") {p <- p + geom_jitter(size=size, height = 0)+
      stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.5, color = "white", size = 1, linetype = "solid")}
    if (plots[i] == "vlnplot") {p <- p + geom_violin()}
  }
  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype="dashed")}
  #Add titles
  if (is.null(main)){
    if(is.meta(var, object)){main <- var}
    if(is.gene(var, object)){main <- paste0("Expression of Wk19 ", var)}
  }
  #Set colors to MYcolors, set the x-axis formatting, and add titles.
  p <- p + scale_fill_manual(values=MYcolors[colors])+ scale_color_manual(values=MYcolors[c(8,2)]) +
    theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2, size=12))+
    ggtitle(main, sub) + xlab(xlab)
  #Change the x-axis labels if requested.
  if (!is.null(labels)) {p <- p + scale_x_discrete(labels=labels)}
  return(p)
}

is.meta <- function(test, object=DEFAULT){test %in% metas(object)}

is.gene <- function(test, object=DEFAULT){
  if(typeof(object)=="character")
  {
    if (classof(object)=="seurat"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    } 
    if (classof(object)=="RNAseq"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    } 
  } else {
    if (class(object)=="seurat"){
      return(test %in% rownames(object@raw.data))
    }
    if (class(object)=="RNAseq"){
      return(test %in% rownames(object@counts))
    }
  }
}

metas <- function(object=DEFAULT){if(typeof(object)=="character"){
  names(eval(expr = parse(text = paste0(object,"@meta.data"))))
} else {names(object@meta.data)}
}

meta <- function(meta = "age", object=DEFAULT){
  if(meta=="ident"){as.character(eval(expr = parse(text = paste0(object,"@ident"))))}
  else{eval(expr = parse(text = paste0(object,"@meta.data$",meta)))}
}

gene <- function(gene, object=DEFAULT){
  if (classof(object)=="seurat"){
    ind <- grep(paste0("^",gene,"$"),rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    OUT <- eval(expr = parse(text = paste0(object,"@data[gene,]")))
  }
  if (classof(object)=="RNAseq"){
    ind <- grep(paste0("^",gene,"$"),rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    OUT <- eval(expr = parse(text = paste0(object,"@rlog[ind,]")))
  }
  OUT
}

extDim <- function(reduction.use, dim=1, object=DEFAULT){
  if (classof(object)=="seurat"){
    OUT <- eval(expr = parse(text = paste0(object,"@dr$",reduction.use,"@cell.embeddings[,",dim,"]")))
  }
  if (classof(object)=="RNAseq"){
    OUT <- eval(expr = parse(text = paste0(object,"@pca[[1]]$x[,",dim,"]")))
  }
  OUT
}

classof <- function (object = DEFAULT){
  class(eval(expr = parse(text = paste0(object)))) 
}

rankedBarcodes <- function(object){
  #Takes in a Seurat object and spits out a cellranger-like rankedBarcodes plot for all the cells in the object
  UMI.counts <- object@meta.data$nUMI[order(object@meta.data$nUMI, decreasing=TRUE)]
  Barcodes <- 1:length(UMI.counts)
  plot(log10(Barcodes),log10(UMI.counts), type="s")
}

fisher.compar <- function(genes, universe, db, enriched = T, p=0.01){
  #This function will run a Fishers Exact Test on the overlap of genes & each of the annotations in db to see
  # if there is enrichment within genes.
  #   genes = a list of gene names.
  #   db = a msigdbr database dataframe.
  #   universe = the universe of genes that genes should be said to come from for enrichment calculation.
  #   enriched = logical, default=TRUE, whether to trim the output to only significant genesets
  #   p = the p_adj cutoff for trimming of significant genesets
  
  #Pull the list names
  lists <- unique(as.character(db$gs_name))
  
  #Pull the genes in each annotation, QUICKLY, as in without a grep.
  #ORIGINAL:list.genes <- sapply(lists, function(X) list(grep(X,as.character(db$gs_name))))
  next.ind <- (1:length(db$gs_name))[!(duplicated(db$gs_name))]
  list.genes <- sapply(1:(length(lists)-1), function(X) list((next.ind[X]):(next.ind[X+1]-1)))
  list.genes <- c(list.genes, list((next.ind[length(lists)]):length(db$gs_name)))
  
  #Initialize the rest.of.universe variable
  rest.of.universe <- universe[!(universe %in% genes)]
  
  #This code obtains the current dataset, used many times below:
  #db$human_gene_symbol[unlist(list.genes[X])]
  
  #### "Populate the table" ####
  # Number of genes in the cluster in the geneset
  A <- sapply(1:length(lists),
              function(X) sum(genes %in% db$human_gene_symbol[unlist(list.genes[X])]))
  # Number of genes in the cluster NOT in the geneset
  B <- sapply(1:length(lists),
              function(X) sum(!(genes %in% db$human_gene_symbol[unlist(list.genes[X])])))
  # Number of genes in the rest of the universe in the geneset
  C <- sapply(1:length(lists),
              function(X) sum(rest.of.universe %in% db$human_gene_symbol[unlist(list.genes[X])]))
  # Number of genes in the rest of the universe NOT in the geneset
  D <- sapply(1:length(lists),
              function(X) sum(!(rest.of.universe %in% db$human_gene_symbol[unlist(list.genes[X])])))
  #### Run the testing ####
  tests <- sapply(1:length(lists),
                  function(X) fisher.test(matrix(c(A[X],B[X],C[X],D[X]), nrow=2, ncol=2))$p.value)
  
  #Make the dataframe
  OUT <- data.frame(geneset = lists, pval = tests, A_overlap=A, expected=((A+C)/(B+D))*(A+B), B=B, C=C, D=D, geneset.size= sapply(1:length(lists), function(X) length(unlist(list.genes[X]))), geneset.in.universe=A+C, stringsAsFactors = F)
  
  #Trim to enriched only (unless input enriched = F)
  if (enriched){ OUT <- OUT[OUT$A_overlap>OUT$expected,]}
  
  #Trim to pval < cutoff.  default = 0.01
  OUT <- OUT[OUT$pval<=p,]
  
  #Order by pvalue
  OUT <- OUT[order(OUT$pval),]
  
  #RETURN
  OUT
}

fisher.mult <- function(genes, universe, dbs, enriched = T, p=0.01){
  #This function will run fisher.compar on a list of name of databases.  The dbs list will also be used
  # to name to dataframes spit out.
  #   genes = a list of gene names to check for enrichment.
  #   dbs = a list of names, in string form, of msigdbr databases already loaded in the memory.
  #   universe = the universe of genes that genes should be said to come from for enrichment calculation.
  #   enriched = logical, default=TRUE, whether to trim the output to only significant genesets
  #   p = the p_adj cutoff for trimming of significant genesets
  OUT <- list()
  OUT <- sapply(1:length(dbs), function(X) list(fisher.compar(genes, universe, db = eval(expr = parse(text = dbs[X])), enriched=enriched, p=p)))
  names(OUT) <- dbs
  OUT
}