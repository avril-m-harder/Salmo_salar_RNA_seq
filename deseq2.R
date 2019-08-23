setwd("/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/stringtie_featureCounts_deseq2_DGE_analyses")
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# biocLite("apeglm")

library("apeglm");library("DESeq2");library("ggplot2");library("dendextend");library("genefilter");library("gplots");library("RColorBrewer");library("rtracklayer");library("beyonce")

#--<>--<>--<>--<>-- Data loading and prep --<>--<>--<>--<>--
## read in count data generated with featureCount
gene.counts <- read.table("stringtie_featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

## read in CSV with sample information
col.data <- read.csv("sample_info.csv", row.names = 1)
## set family ID as a factor
col.data[,2] <- as.factor(col.data[,2])

## check to make sure that sample names are in the same order in the gene count table and the
## sample info table
all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>

################### ANALYSIS ##########################################
#--<>--<>--<>--<>-- Construction of DESeqDataSet object and DE analysis --<>--<>--<>--<>--
# test for effects of treatment while controlling for the effect of family; dds = DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ family + treatment)

dds$treatment <- relevel(dds$treatment, "untreated") ## set the untreated group as the control
dds

## renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG wherever gene name is known)
assembly <- readGFF("stringtie_all_merged.annotated.gtf")             ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)

## create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
gene.name <- assembly$gene_name[gene_idx]
transcript.id <- assembly$transcript_id[gene_idx]
xloc <- assembly$xloc[gene_idx]
gene_names <- cbind(gene.name, transcript.id, xloc)

## adds unique gene names to dds
dds@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep="|")
## check to make sure that no gene names in dds are duplicates
which(duplicated(dds@rowRanges@partitioning@NAMES))

## DESeq() = differential expression analysis based on the negative binomial distribution
dds <- DESeq(dds)

## examine distribution of dispersion values
plotDispEsts(dds, ylim=c(1e-6, 1e1))
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--

#--<>--<>--<>--<>-- Extract and examine results from DESeq analysis --<>--<>--<>--<>--
res <- results(dds, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)

## examine MA plot (normalized counts vs. log2fold changes)
plotMA(res, ylim=c(-3,3), xlim=c(1e-2,1e6))

## examine histograms of p-values and adjusted p-values
hist(res$pvalue, breaks=20, col="grey")
hist(res$padj, breaks=20, col="grey")
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

################### DATA MANIPULATION AND EXPLORATION ############################
#--<>--<>--<>--<>-- Log transformation and distance calculation --<>--<>--<>--<>--
## rlog() = transforms count data to the log2 scale in a way that minimizes differences between
## samples for genes with small counts, and which normalizes with respect to library size
rld <- rlog(dds)
head(assay(rld))
save.image("through_logxform.RData")

## calculate Euclidean distances between all samples to examine overall similarity
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

## construct heatmaps to examine sample patterns with different clustering methods
rownames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
colnames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
dists.ward <- hclust(sampleDists, method="ward.D2")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.ward),
          margins=c(6,7))
dists.avg <- hclust(sampleDists, method="average")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.avg),
          margins=c(6,7))
heatmap.2(sampleDistMatrix, trace="none", dendrogram="column",
          scale="row", margins=c(6,7))

#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

## for plotPCA(), >getMethod("plotPCA","DESeqTransform") to see source code

# --<>--<>--<>--<>-- PCA to examine effect of treatment --<>--<>--<>--<>--<>--<>--
plot.treat.data <- plotPCA(rld, intgroup = c("treatment"), returnData=TRUE, n=500) ## uses n most variable genes -- not necessarily sig DEGs!
percentVar <- round(100*attr(plot.treat.data, "percentVar"))
my.colors <- beyonce_palette(94,2,type=c("discrete"))
ggplot(plot.treat.data, aes(PC1, PC2, color=treatment)) +
  scale_color_manual(values=c(my.colors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
# --<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

#--<>--<>--<>--<>-- PCA to examine effects of treatment + family --<>--<>--<>--<>--
plot.all.data <- plotPCA(rld, intgroup = c("family", "treatment"), returnData=TRUE, n=500) ## uses n most variable genes -- not necessarily sig DEGs!
percentVar <- round(100*attr(plot.all.data, "percentVar"))
my.colors <- beyonce_palette(66,9,type=c("continuous"))
ggplot(plot.all.data, aes(PC1, PC2, color=family, shape=treatment)) +
  scale_color_manual(values=c(my.colors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--

################### ID OF DEGs ############################################
#--<>--<>--<>--<>-- How many and which genes are sig DE? --<>--<>--<>--<>--
sum(res$padj < 0.05, na.rm=T)            
sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm=T)
resSig <- res[which(res$padj < 0.05),]   ## put all genes w/ padj < 0.05 in resSig
plotMA(resSig, ylim=c(-3,3), xlim=c(1e-2,1e6)) ## examine MA plot 
head(resSig[order(resSig$log2FoldChange),])    ## sort by log2foldchange and show most down-regulated
tail(resSig[order(resSig$log2FoldChange),])    ## sort by log2foldchange and show most up-regulated
write.csv(as.data.frame(resSig), file="p_05_degs.csv")     ## write CSV with all DEGs, padj < 0.05
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--

################### VISUALIZING DEGs #######################################
#--<>--<>--<>--<>-- Heatmap w/ top significant DEGs --<>--<>--<>--<>--

topSigGenes <- head(resSig[order(resSig$padj),],n=10)      ## get 10 genes with lowest padj values
colours <- beyonce_palette(64,300,type=c("continuous"))

#!# genes on y, samples on x
heatmap.2(assay(rld)[rownames(topSigGenes),], 
          scale="row", ## plot heatmap using these genes
          trace="none", dendrogram="column",margins=c(4,10), key.title=NA,
          cexRow=1, cexCol=1,
          col = colours)
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

#--<>--<>--<>--<>-- PCA plots - DEGs by cut-offs --<>--<>--<>--<>--
p.cutoff <- 0.05
fc.cutoff <- 0

topSigGenes <- res[which(res$padj < p.cutoff & abs(res$log2FoldChange) >= fc.cutoff),]
num.genes <- topSigGenes@nrows
rownames(rld)
rownames(topSigGenes)
rld.Sig <- rld[rownames(rld) %in% rownames(topSigGenes)]

plot.all.data <- plotPCA(rld.Sig, intgroup = c("family", "treatment"), returnData=TRUE, ntop=num.genes)
percentVar <- round(100*attr(plot.all.data, "percentVar"))
my.colors <- beyonce_palette(66,9,type=c("continuous"))

## color PCA by family
ggplot(plot.all.data, aes(PC1, PC2, color=family, shape=treatment)) +
  ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +       
  scale_color_manual(values=c(my.colors)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

## color PCA by treatment
ggplot(plot.all.data, aes(PC1, PC2, color=treatment, shape=treatment)) +
  ggtitle(paste("Top",num.genes,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +       
  scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
#--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--
