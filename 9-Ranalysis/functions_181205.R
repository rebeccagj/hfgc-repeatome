#These are my (sc)RNAseq functions and some things that they rely on.
MYcolors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
library(ggplot2, quietly = T)
library(dplyr, quietly = T)

# ### Bulk RNASeq block ####
# library(DESeq2, quietly = T)
# 
# Class <- setClass("RNAseq",
#                   representation(
#                     counts = "matrix",
#                     dds = "DESeqDataSet",
#                     rlog = "matrix",
#                     meta.data = "data.frame",
#                     pca = "list",
#                     var.genes = "character",
#                     samples = "character",
#                     exp.filter = "logical",
#                     CVs = "numeric"
#                     ),
#                   prototype(
#                     counts = matrix(),
#                     dds = new("DESeqDataSet"),
#                     rlog = matrix(),
#                     meta.data = data.frame(),
#                     pca = list(),
#                     var.genes = character(),
#                     samples = character(),
#                     exp.filter = logical(),
#                     CVs = double()
#                   )
#   )
# 
# #setMethod("initialize", "RNAseq", 
# #          function(.Object, ...) {
# #            .Object <- callNextMethod()
# #            if(length(.Object@x) != length(.Object@y))
# #              stop("specified x and y of different lengths")
# #            .Object
# #          })
# 
# NewRNASeq <- function(dds, #A DESeq object, *the output of DESeq()*
#                       populate_PCA = FALSE,#If changed to TRUE, function will:
#                                            # auto-populate the rlog, var.genes, and PCA fields.
#                       Ngenes = 2500, #How many genes to use for running PCA, (and how many genes
#                                      # will be stored in @var.genes)
#                       blind = FALSE #Whether or not the rlog estimation should be blinded to sample info.
#                                     # Run `?rlog` for more info
#                       ){
#   #INPUTS
#   #dds                The DESeq2 object for your data, *the output of DESeq()*
#   #populate_PCA       FALSE by default.  Whether @rlog @var.genes and @pca[[1]] population are desired.
#   #Ngenes             How many genes to use for the PCA
#   #blind              Whether rlog estimation should be blinded to sample info. Run `?rlog` for more info
#   
#   ########## Create the Object ########################
#   #Create the object with whatever inputs were given, a.k.a. creates objects@counts and any other level
#   # within str(object).  Will all be NULL unless provided in the function call
#   object <- new("RNAseq", dds = dds)
#   
#   ########## Run Autopopulations ######################
#   # Will run by default because this function requires a dds object to be given.
#   ##populate dds
#   object@dds <- dds
#   ##populate @counts
#   object@counts <- counts(dds)
#   ##populate @samples
#   object@samples <- colnames(object@counts)
#   ##populate some of @meta.data
#   #   1st add samples, then Nreads.
#   object@meta.data <- data.frame(Samples = object@samples,
#                                  Nreads = colSums(object@counts))
#   
#   ##Also add colData from dds to @meta.data slot
#   #Turn colData into a data.frame, and merge that with current meta.data, BUT do not include any
#   # dublicate sets.  For example, Samples will be ignored in colData because it was already grabbed
#   # from the counts matrix
#   object@meta.data <- cbind(object@meta.data,
#                             data.frame(object@dds@colData@listData)[(!duplicated(
#                               c(names(data.frame(object@dds@colData@listData)),names(object@meta.data)),
#                               fromLast=T
#                               ))[1:length(object@dds@colData@listData)]])
#   ########## Will run if populate_PCA = TRUE ##################
#   ##Will populate: rlog, pca, var.genes
#   if (populate_PCA){
#     ##populate rlog
#     object@rlog <- assay(rlog(object@dds, blind = blind))
#     #Filter rlog to only the genes expressed in at least 75% of samples from test group (ONLY WORKS FOR ONE TEST GROUP)
#       test_meta <- strsplit(as.character(object@dds@design), split = "~")[[2]][1]
#       #Store this metadata as an easily accessible variable to speed up the next step.
#       classes <- meta(test_meta)
#       #For each gene, return TRUE if... the gene is expressed in >75% of samples from each condition used to build the dds
#       ##populate exp.filter
#       object@exp.filter <- sapply(1:dim(object@counts)[1], function(X)
#         #Ensures that the #of classifications = the number TRUEs in what the nested sapply produces
#         length(levels(as.factor(classes)))==sum(
#           #For each classification of the test variable, check for >= 75% expression
#           #This half sets a variable to each of the saparate classifications,
#           # and says to run the next lines on each
#           sapply(levels(as.factor(classes)), function(Y)
#             #This part of the function determine how many of the samples express the gene
#             (sum(object@counts[X, classes==Y]>0))
#             #This part compares the above to (the number of samples of the current classification*75%)
#              >= (sum(classes==Y)*.75)
#           )
#         )
#       )
#     data_for_prcomp <- as.data.frame(object@rlog)[object@exp.filter,]
#     #calculate CV by dividing mean by sd
#     ## populate CVs
#     object@CVs <- apply(X = object@rlog, MARGIN = 1, FUN = sd)/apply(X = object@rlog, MARGIN = 1, FUN = mean)
#     #Trim rlog data and RawCV_rlog by expression filter variable = object@exp.filter
#     #arrange by CV_rank, higher CVs first
#     data_for_prcomp<- data_for_prcomp[order(object@CVs[object@exp.filter], decreasing = T),]
#     ##populate var.genes
#     object@var.genes <- rownames(data_for_prcomp)[1:(min(Ngenes,dim(data_for_prcomp)[1]))]
#     ##populate pca : Run PCA on the top Ngenes CV genes that survive a 75% expression per condition filter
#     object@pca <- list(prcomp(t(data_for_prcomp[1:Ngenes,]), center = T, scale = T))
#   }
#   #OUTPUT: (This is how functions "work" in R.  The final line is what they return.)
#   object
# }

# #####################

#This is Burt lab specific... my labmates generally make a theme that is stored as prettyplot. ####
# I put this in as an example for them.
prettyplot.1 <- theme(text = element_text(size = 14, color="black"),
                      #legend.position="none",
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      axis.text=element_text(color="black"),
                      legend.title = element_text(colour="black", size=14,face="bold"),
                      legend.text = element_text(colour="black", size = 12, face="plain"),
                      plot.background = element_rect(fill = "transparent",colour = NA))


############# MAIN FUNCTIONS ##############
# DBDimPlot ####
DBDimPlot <- function(var="ident", object = DEFAULT, reduction.use = NA, dim.1 = 1, dim.2 = 2, main = NULL,
                      sub = NULL, cells.use = NULL, theme = NA, do.label = F, label.size = 5, highlight.labels = F,
                      show.others=TRUE, size=1, shape=16, low = "#F0E442", high = "#0072B2", colors = 1:8,
                      range = NULL, auto.title = T, ellipse = F, xlab=NA, ylab = NA, legend.size = 5){
  #Makes a plot where color is overlayed onto the dim.reduction plot of choice.
  #
  #object                 the Seurat Object
  #var                    Target Variable = either metadata name or values
  #reduction.use          "pca", "tsne", "ica"
  #dim.1                  If dim.reduction plot, x-axis, component number if reduction.use!=null
  #dim.2                  If dim.reduction plot, y-axis, component number if reduction.use!=null
  #main                   plot title
  #sub                    plot subtitle
  #cells.use              cells to show
  #show.others            TRUE by default, whether other cells should be shown in the background
  #size                   number for size of all highlighted points
  #shape                  number for setting shape OR name of metadata to use for setting shape
  #low                    color for lowest values of var
  #high                   color for highest values of var
  #colors                 indexes of color from MYcolors
  #range                  limits for the color scaling
  #theme                  Allows setting of a theme. Uses theme_bw if not provided.
                          # To provide, either say "prettyplot" and the prettyplot.1 declared in
                          # this file will be used.  Or provide as full "list" form.
  #xlab & ylab            labels for the x and y axes.
  #do.label               Whether to add text labels at the center (median) of clusters for grouping vars
  #label.size             Size of the the labels text
  #legend.size            The Size to increase the plotting of legend shapes to 
  
  #Establish Defaults
  #If cells.use = NA (was not provided), populate it to be all cells or all samples.
  if (classof(object)=="seurat" & is.null(cells.use)) {cells.use <- eval(expr = parse(text = paste0(object,"@cell.names")))}
  if (classof(object)=="RNAseq" & is.null(cells.use)) {cells.use <- meta("Samples", object = object)}
  #If reduction.use = NA (was not provided), populate it to be tsne or pca.
  if (classof(object)=="seurat" & is.na(reduction.use)) {reduction.use <- "tsne"}
  if (classof(object)=="RNAseq" & is.na(reduction.use)) {reduction.use <- "pca"}

  #Build data for populating dat, the data.frame for plotyting.
  #Determine the identity of the provided 'var' and populate Y, the variable used for coloring.
  if(typeof(var)=="character"){
    #If "ident" pull the @ident object from the seurat object
    if(var == "ident"){Y <- eval(expr = parse(text = paste0(object, "@ident")))}
    #If "is.meta" pull the @meta.data$"var" from the RNAseq or seurat object
    if(is.meta(var, object)){Y <- meta(var, object)}
    #If "is.gene" pull the gene expression data from the RNAseq or seurat object
    if(is.gene(var, object)){Y <- gene(var, object)}
    #Otherwise, var is likely a full set of data already, so just make Y = var
  } else {Y <- var}
  #Determine the identity of the provided 'shape'.
    #If it is a meta.data name, pull the meta.  else (it is a number) just carry it through
  if(typeof(shape)=="character"){Shape <- meta(shape, object)} else{Shape <- shape}
  
  #Populate the data.frame to be used for plotting
  full_dat <- data.frame(Y = Y,
                         dim1 = extDim(reduction.use,dim.1,object),
                         dim2 = extDim(reduction.use,dim.2,object),
                         size = size,
                         shape = Shape,
                         Group = Y)

  #Subset to cells.use
  if (classof(object)=="seurat"){
    Others_dat <- full_dat[!(eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use),]
    Target_dat <- full_dat[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use,]
  } 
  if (classof(object)=="RNAseq"){
    Others_dat <- full_dat[!(meta("Samples") %in% cells.use),]
    Target_dat <- full_dat[meta("Samples") %in% cells.use,]
  } 
  
  ###Start building the plot###
  p <- ggplot() +
    #Remove the legend title
    theme(legend.title=element_blank())
    #If not provided, add default axis labels (ex. "tsne1" or "pca2")
    if (is.na(xlab)) {p <- p + xlab(paste0(reduction.use,dim.1))
    } else {p <- p + xlab(xlab)}
    if (is.na(ylab)) {p <- p + ylab(paste0(reduction.use,dim.2))
    } else {p <- p + ylab(ylab)}
  #Then Add more layers:
  
  ###Add the data###
  #Make gray dots on the bottom layer if show.others = T and cells.use is a subset of all the cells / samples.
  if (show.others & dim(Others_dat)[1]>1) {
    p <- p + geom_point(data=Others_dat, aes(x = dim1, y = dim2), size=0.5, color = "gray90")
  }
  #Overlay the target data on top
  #If 'shape' input was the name of a meta.data, aka type=character, treat shape as an aesthetic for performing grouping.
  # Otherwise it is a number and belongs outside of aes.  
  if (typeof(shape)=="character") {
    p <- p + geom_point(data=Target_dat, aes(x = dim1, y = dim2, colour = Y, shape= shape), size=size)  
  }  else {
    p <- p + geom_point(data=Target_dat, aes(x = dim1, y = dim2, colour = Y), shape= shape, size=size)
  }
  
  ###Add and ellipse###
  ### Draw an ellipse if ellipse = T.
  if (ellipse) { p <- p + stat_ellipse(data=Target_dat,
                                       aes(x = dim1, y = dim2, colour = Y),
                                       type = "t",
                                       linetype = 2,
                                       size = 0.5
                                       )}
  
  ###Add titles###
  #If not provided, autogenerate based on the identity of var
  if (is.null(main) & auto.title==T){
    #If var is a meta.data, make the title = var
    if(is.meta(var, object)){main <- var}
    #If var is a gene, make the title "Expression of "var
    if(is.gene(var, object)){main <- paste0("Expression of ", var)}
  }
  #If main was provided, use that as the main title
  if (!is.null(main)) {
    p <- p + ggtitle(main, subtitle = sub)
  }
  #If sub was provided, use that as the subtitle
  # if (!is.null(sub)) {
  #   p <- p + ggtitle(subtitle = sub)
  # }
  
  ### Add Labels ###
  if (do.label) {
    #Make a text plot at the median x and y values for each cluster
    #Determine medians
    cent.1 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim1[Target_dat$Y==X]))
    cent.2 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim2[Target_dat$Y==X]))
    #Add labels
    if (highlight.labels){
      #Add labels with a white background
      p <- p + 
        geom_label(data = data.frame(x=cent.1, y=cent.2),
                  aes(x = x, y = y, label = levels(as.factor(Y))),
                  size = label.size)
    } else {
      #Add labels without a white background
      p <- p + 
        geom_text(data = data.frame(x=cent.1, y=cent.2),
                  aes(x = x, y = y, label = levels(as.factor(Y))),
                  size = label.size)
    }
  }
  
  ### Set the colors ###
  ### Also change the size of the dots in the legend if showing groupings ###
  #If var yielded a list of groups for plotting (should be in the form of a list of strings = character, or a factor = integer)
  if (typeof(Y)=="character" | typeof(Y)=="integer"){
    #If the number of levels/groups is less than nine, use my colors set.
    if (length(levels(as.factor(as.character(Target_dat$Y))))<9){
      p <- p + scale_colour_manual(values = MYcolors[colors])
      #Also change the size of the dots in the legend unless legend.size was set to NA
      if (!is.na(legend.size)){
        p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
      }
    }
  } else {
    #Otherwise, the data is continous, so set a gradient that goes from 'low' input color to 'high' input color.
    p <- p + scale_color_gradient(low= low, high = high, limits = range)
  }
  
  ### Set the theme ###
  #Use theme_bw if theme = NA (was not provided), use prettyplot.1 if "prettyplot" was provided, or
  # use provided theme if a full one is provided, aka = a list.  *Does not check if your theme is complete.
  if (is.na(theme)){
    p <- p + theme_bw() + theme(legend.title=element_blank())
  } else {
    if (theme=="prettyplot") {p <- p + prettyplot.1}
    if (typeof(theme)=="list") {p <- p + theme}
  }
  
  #DONE, return the plot as output
  return(p)
}

# DBPlot ####
DBPlot <- function(var, object = DEFAULT, main = NULL, sub = NULL, group.by = "Tage",
                   color.by = "age", cells.use = NULL,
                   size=0.1, shape=16, low = "#F0E442", high = "#0072B2", colors = c(1:8),
                   ylab = NULL, xlab = NULL, plots = c("jitter","vlnplot"), hline=NULL, labels = NULL,
                   range = NULL, jitter.width=0.2, boxplot.width = 0.2, jitter.color = "black",
                   boxplot.color = "black", show.outliers = F, color.panel = MYcolors, y.breaks = NA){
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
  #jitter.width           the width/spread of the jitter in the x direction
  #boxplt.width           the width/spread of the boxplot in the x direction
  #color.panel            the set of colors to draw from
  
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
  p <- ggplot(Target_dat, aes(x=sample, y=Y, fill=color)) +
    #And set the legend to not have a title
    theme(legend.title=element_blank())
  #Add the y label
  p<- p + ylab(ylab)
  #Set the y-axis limits if a range is given.
  if (!is.null(range)){p <- p + ylim(range)}
  ###Add data based on what is requested in plots, ordered by their order
  for (i in 1:length(plots)){
    if (plots[i] == "boxplot") {p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                                                      outlier.shape = ifelse(show.outliers,16,NA))} 
    if (plots[i] == "jitter") {p <- p + geom_jitter(size=size, width=jitter.width, height = 0, color = jitter.color)}
    if (plots[i] == "vlnplot") {p <- p + geom_violin()}
  }
  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype="dashed")}
  #Add titles
  if (is.null(main)){
    if(is.meta(var, object)){main <- var}
    if(is.gene(var, object)){main <- paste0("Expression of ", var)}
  }
  #Set colors to MYcolors, set the x-axis formatting, and add titles.
  p <- p + scale_fill_manual(values=MYcolors[colors])+
    #scale_color_manual(values=MYcolors[colors]) +
    theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2, size=12))+
    ggtitle(main, sub) + xlab(xlab)
  #Change the x-axis labels if requested.
  if (!is.null(labels)) {p <- p + scale_x_discrete(labels=labels)}
  #Change y-axis limits/breaks if requested
  if (!is.na(y.breaks)) {
    p <- p + scale_y_continuous(breaks = y.breaks) + coord_cartesian(ylim=c(min(y.breaks),max(y.breaks)))
  }
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