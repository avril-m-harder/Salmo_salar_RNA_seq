### Weighted gene co-expression network analysis (WGCNA) for RNA-seq data ###
### Largely based on WGCNA tutorials at: 
### https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
setwd("~/Documents/rna_seq/seq_data_processing_notes_and_analyses/wgcna/")

library("WGCNA")

##### DATA LOADING AND PREP #####
# load("through_logxform.RData") ## saved during DESeq2 analysis

new.dds <- estimateSizeFactors(dds)
counts <- counts(new.dds, normalized=TRUE)

options(stringsAsFactors = FALSE)
dim(counts) ## 36 samples, 94,785 transcripts
rownames(counts) ## sample names

## check for genes and samples with too many missing values
gsg <- goodSamplesGenes(counts, verbose = 3)
gsg$allOK ## <-- this should return "TRUE" ==> if it doesn't, run below loop

# if (!gsg$allOK)
# {
#   ## this is where genes that weren't "good" were filtered out.
#   ## optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(counts)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(counts)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   counts = counts[gsg$goodSamples, gsg$goodGenes]
# }

## how many genes have at least 30 samples with non-zero counts?
dim(counts[rowSums(counts==0) <= 30,])
## make this the new criteria for counts to speed things up for outlier check step
counts <- counts[rowSums(counts==0) <= 30,]

# ## cluster samples to see if there are any obvious outliers
# sampleTree = hclust(dist(counts), method = "average")
# par(cex = 0.6)
# par(mar = c(2,4.5,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# ## Static tree cutting
# # Plot a line to show the cut
# abline(h = 80000, col = "red", lty=3);
# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
# table(clust) ## how many will be cut [0] and how many will be kept [1]
# # clust 1 contains the 35 samples we want to keep.
# keepSamples = (clust==1)
# # datExpr = counts[keepSamples, ]
datExpr <- t(counts)
# sampleTree = hclust(dist(datExpr), method = "average")
# par(cex = 0.6)
# par(mar = c(2,4.5,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
## gene expression data now ready for expression analysis

## read in trait data and match info. to expression samples
## stack and treatment ID recoded as #s
traitData <- read.csv("refined_fry_traits.csv")
dim(traitData)
names(traitData)

## make data frame to match expression data
sampnames <- rownames(datExpr)
traitrows <- match(sampnames, traitData$fry.id)
datTraits <- traitData[traitrows, -c(1,2)]
rownames(datTraits) <- traitData[traitrows, 1]
collectGarbage();

## important data are now in datExpr and datTraits
## want to see how the traits relate to the sample dendrogram
sampleTree2 <- hclust(dist(datExpr), method="average")
## convert traits to a colr representation; white = low values, red = high value, grey = missing value
traitColors <- numbers2colors(datTraits, signed=FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=names(datTraits),
                    main="Sample dendrogram and trait heatmap")

##### Network Construction #####
# choose a set of soft-thresholding powers to start with
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
## call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## plot the results
par(mfrow = c(1,2))
cex1 = 0.9
## scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
## going with 6 

## construct gene network and ID modules
## power determined above
## TOM = topological overlap matrix
## blockwiseModules() == automatic module detection via dynamic tree cutting
net <- blockwiseModules(datExpr, power = 6,
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.10,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        networkType = "signed",
                        saveTOMFileBase = "rep4_TOM", 
                        verbose = 3)
## check to see how many modules were ID'ed and the size of the modules
## IF THIS ERRORS OUT: run search(); ".GlobalEnv" and "package:WGCNA" should be the first 2 in the list
## if they're not, restart R and only load WGCNA

# save.image("rep4_allRData.RData")

##### STARTING POINT FOR AUTO MODULE ID W/ POWER=6, MINMODSIZE=30 , MERGECUTHEIGHT=0.10 #####

# load("rep4_allRData.RData")

table(net$colors)

## convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
## plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## save the module assignment and module eigengene information necessary for subsequent analysis
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];

##### Relate modules to external info and ID important genes #####

## define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

## recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# var.explain <- moduleEigengenes(datExpr, moduleColors)$varExplained
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

##### GETTING LIST OF GENE IDS #####
sig.modules <- moduleTraitPvalue[apply(moduleTraitPvalue, 1, function(row) {any(row < 0.05)}),]
list.sig.modules <- rownames(sig.modules)
for (i in 1:length(list.sig.modules)) {
  list.sig.modules[i] <- substring(list.sig.modules[i], 3) ## remove "ME" from the beginnings of all module names
}

module.gene.IDs <- data.frame()
i <- 1
for (module in list.sig.modules) {
  tmp.gene.IDs <- data.frame()
  mod.genes <- dimnames(datExpr)[[2]][moduleColors==module]
  tmp.gene.IDs[i:(length(mod.genes)),1] <- module
  tmp.gene.IDs[i:(length(mod.genes)),2] <- mod.genes
  module.gene.IDs <- rbind(module.gene.IDs, tmp.gene.IDs)
}
module.gene.IDs
sum(table(module.gene.IDs))

# write.csv(module.gene.IDs, "/Users/Avril/Desktop/rep4_module_gene_IDs.csv")

##### BUILD LIST OF MODULES SIG. ASSOC. WITH VARS, GET CORS AND P-VALUES #####
head(moduleTraitCor)
head(moduleTraitPvalue)

all(rownames(moduleTraitCor) == rownames(moduleTraitPvalue))
all(colnames(moduleTraitCor) == colnames(moduleTraitPvalue))

p.vals <- data.frame(module=character(),
                     n.genes=numeric(),
                     pred.var=character(),
                     cor=numeric(),
                     p.value=numeric(),
                     stringsAsFactors=FALSE)
temp <- data.frame(module=character(),
                   n.genes=numeric(),
                   pred.var=character(),
                   cor=numeric(),
                   p.value=numeric(),
                   stringsAsFactors=FALSE)

for(row in 1:nrow(moduleTraitPvalue)) {
  for (col in 1:ncol(moduleTraitPvalue)) {
    if (moduleTraitPvalue[row,col]<0.05) {
      temp.mod <- substring(rownames(moduleTraitPvalue)[row],3)
      temp[1,1] <- rownames(moduleTraitPvalue)[row]
      temp[1,3] <- colnames(moduleTraitPvalue)[col]
      temp[1,4] <- moduleTraitCor[row,col]
      temp[1,5] <- moduleTraitPvalue[row,col]
      temp[1,2] <- length(colnames(datExpr)[moduleColors==temp.mod])
      p.vals <- rbind(p.vals, temp)
    }
  }
}

# write.csv(p.vals,"/Users/Avril/Desktop/rep4_sig_modules_cors_pvals.csv", row.names=F)

##### GRAB "TOP" SIG. MODULE:VARIABLE ASSOCIATIONS #####
bonf <- (0.05/length(rownames(moduleTraitCor)))
top.modules <- p.vals[which(p.vals$p.value <= bonf & p.vals$pred.var == "treatment"),]
rownames(top.modules) <- top.modules[,1]
top.modules <- top.modules[,-1]