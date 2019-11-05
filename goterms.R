## Script for running GO term enrichment analysis on a set of genes using a reference
## set of GO terms (Ssal_ICSASG_v2_GOAccession.txt) published with the reference
## Atlantic salmon genome
##
## category_gene_go_info.csv includes column of gene names, column of corresponding
## GO terms, and column denoting which group of genes each gene belongs to
##
## code written by J. Willoughby

# source("http://bioconductor.org/biocLite.R")
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")
# mc <- "https://cran.r-project.org/src/contrib/Archive/metacoder/metacoder_0.1.3.tar.gz"
# install.packages(mc, repos=NULL, type="source")

library(topGO)
library(Rgraphviz)

setwd("~/analyses_ssalar_rna_seq/category_go_analysis")

#read in genes of interest
interest_goALL = read.table("category_gene_go_info.csv", header=T, sep=",")
colnames(interest_goALL)[2] <- "module"
toselect = function(data){return(names(data[data==1]))} #selection function

#read in genome-wide information and set up objects for further processing
genome_go = read.table("Ssal_ICSASG_v2_GOAccession.txt", sep="\t", header=TRUE, quote="", as.is=TRUE, comment.char="")
genome_annots = cbind(genome_go$SeqName, genome_go$GO.ID)
genome_go = cbind(genome_go$SeqName, genome_go$GO.ID, rep(0, nrow(genome_go)))
colnames(genome_go) = c("seqID", "GO", "select")

#prepare gene-to-GO annotations - takes forever so can read in file if available
# ugenes = unique(genome_annots[,1])
# ugenes_annots = matrix(nrow=length(ugenes), ncol=2)
# for(g in 1:length(ugenes)){
#   gotemp = NULL
#   temp = genome_annots[genome_annots[,1]==as.character(ugenes[g]),,drop=FALSE]
#   if(nrow(temp)>1){
#     gotemp = paste(t(temp[,2]), collapse=",")
#   }else if(nrow(temp)==1){
#     gotemp = temp[1,2]
#   }
#   ugenes_annots[g,1] = temp[1,1]
#   ugenes_annots[g,2] = gotemp
# }
# ugenes_annots = ugenes_annots[1:g,]
# write.table(ugenes_annots, "ugenes_annots.txt", sep="\t", col.names=F, row.names=F)
backup_ugenes_annots = as.matrix(read.table("ugenes_annots.txt", sep="\t", header=F))

#iterate over modules
mods = unique(interest_goALL$module)
sumGO = data.frame(module = mods, elim = rep(NA, length(mods)), weight01 = rep(NA, length(mods)), elimweight01 = rep(NA, length(mods)))
for(g in 1:length(mods)){
  ugenes_annots   = backup_ugenes_annots
  interest_go     = interest_goALL[interest_goALL$module==as.character(mods[g]),]
  interest_go     = interest_go[,2:ncol(interest_go),]
  interest_go_out = NULL
  for(i in 1:length(interest_go[,3])){
    gotemp = NULL
    temp = interest_go[interest_go[,3]==as.character(interest_go[i,3]),,drop=FALSE]
    if(nrow(temp)>1){
      gotemp = paste(t(temp[,4]), collapse=", ")
    }else if(nrow(temp)==1){
      gotemp = temp[1,2]
    }
    interest_go_out = rbind(interest_go_out, c(as.character(temp[1,3]), gotemp))
  }
  go = cbind(as.vector(interest_go_out[,1]), as.character(interest_go_out[,2]), rep(1, nrow(interest_go_out)))
  go = go[go[,2]!="",,drop=FALSE]
  colnames(go) = c("seqID", "GO", "select")
  
  #save data into  one list
  data = rbind(go, genome_go)
  data = as.data.frame(data)
  
  #pull out gene data info, including names
  genes = as.factor(data[,3])
  names(genes) = data[,1]
  
  #finish up with annotation file
  ugenes_annots = cbind(ugenes_annots, rep(0, nrow(ugenes_annots)))
  allannot = rbind(go, ugenes_annots)
  write.table(allannot, "allannot.map", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  #set up needed data structures
  geneID2GO = readMappings(file = "allannot.map")
  
  #analysis begin
  GOdata.BP = new("topGOdata", description = "Ssalar", ontology = "BP", allGenes = genes, geneSel = toselect, nodeSize = 5, annotationFun=annFUN.gene2GO, gene2GO=geneID2GO)
  
  #run test
  resultF.elim.BP     = runTest(GOdata.BP, algorithm = "elim",     statistic = "fisher") #
  resultF.weight01.BP = runTest(GOdata.BP, algorithm = "weight01", statistic = "fisher") #
  
  #calculate BH scores
  BH.elimF     = data.frame(GO.ID = names(p.adjust(score(resultF.weight01.BP), method = "BH")), Felim.BH     = p.adjust(score(resultF.elim.BP),     method = "BH"))
  BH.weight01F = data.frame(GO.ID = names(p.adjust(score(resultF.weight01.BP), method = "BH")), Fweight01.BH = p.adjust(score(resultF.weight01.BP), method = "BH"))
  BH = merge(BH.elimF, BH.weight01F, by="GO.ID", all.x=TRUE, all.y=TRUE )
  
  #compare results from methods
  top.results = GenTable(GOdata.BP, weight01F = resultF.weight01.BP, elimF  = resultF.elim.BP, orderBy = "elimF", ranksOf = "weight01F", topNodes = 5000)
  
  #write out
  top.results = cbind(top.results, module=rep(mods[g], nrow(top.results)))
  if(g==1){write.table(top.results, "topGOenrichmentresultsBP.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
  }else{write.table(top.results, "topGOenrichmentresultsBP.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)}
  remove(BH.elimF, BH.weight01F, BH)
  
  #find terms that have small p-values for weight01 and elim 
  top.results$weight01F = as.numeric(as.character(top.results$weight01F))
  top.results$weight01F[is.na(top.results$weight01F)] = 0
  top.results$elimF = as.numeric(as.character(top.results$elimF))
  top.results$elimF[is.na(top.results$elimF)] = 0
  topterms = top.results[top.results$weight01F <= 0.001 & top.results$elimF <= 0.001, 1]
  
  #summarize
  sumGO$elim[g]         = length(top.results$elimF[top.results$elimF <= 0.001])
  sumGO$weight01[g]     = length(top.results$weight01[top.results$weight01 <= 0.001])
  sumGO$elimweight01[g] = length(topterms)
  
  #find genes with top go terms
  tgo_genes = NULL
  if(length(topterms)>0){
    for(tt in 1:length(topterms)){
      tgo  = topterms[tt]
      temp = genome_go[genome_go[,2]==as.character(tgo),1,drop=FALSE]
      tgo_genes = temp
      tgo_genes = unique(tgo_genes)
      
      #add module info for later
      tgo_df = data.frame(gene = tgo_genes, GOterm = rep(topterms[tt], length(tgo_genes)), module = rep(mods[g], length(tgo_genes)))
      
      #write out results
      if(g==1){write.table(tgo_df, "topgoterms_genesBP.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
      }else{write.table(tgo_df,    "topgoterms_genesBP.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)}
    }
  }
  
  #generate some stupid and rough plots
  pdf(paste("BPtopterms", mods[g], ".pdf", sep=""), width=20, height=10, useDingbats=FALSE, onefile=TRUE)
  try(showSigOfNodes(GOdata.BP, score(resultF.elim.BP)[names(score(resultF.elim.BP)) %in% topterms], firstSigNodes = length(topterms), useInfo="def"), silent=TRUE)
  dev.off()
  
}
write.table(sumGO, "summarySigGO.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
