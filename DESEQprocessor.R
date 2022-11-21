#! /usr/bin/Rscript
#YuEnoch 13-11-2022
#DESEQprocessor.R 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/
#                and Dr. Matt Holding

#Purpose: uses Merged Counts across all treatments for DESeq2 analysis to identify significantly enriched peptides
#         uses FDR p<0.05 (default) and Bonferroni p<0.05 to be defined as significant

#         Other Plots:
#         Plots most sigificant gene in each treatment
#         Plots comparison in counts between FDR and Bonferroni

#Changes across Experiments: alter accordingly
# 1. Mutagenesis for 5 Amino Acids
# 2. NNK Mutagenesis

#Edit the Following: working directory, treatment columns, replicates, Merged DESeq Csv file location
#change folder to the location of your folder
folder = "C:/Users/enoch/Downloads/Github ALL"

coldata = c("cathepsinG","cathepsinG","hpr3","hpr3","elastase","elastase","N","N")
replicates = c("C1","C2","H1","H2","E1","E2","N1","N2")
treatments <- c("cathepsinG", "hpr3", "elastase")
experiment = "Neutrophils"
data <- read.csv(paste0(folder, "/", experiment, "_merged.csv"),header=TRUE,row.names=1)

setwd(folder)


#coldata = c("control","control","control","AEBSF","AEBSF","AEBSF","EDTA","EDTA","EDTA","AEBSF+EDTA","AEBSF+EDTA","AEBSF+EDTA", "N","N","N")
#replicates = c("C1",    "C2",   "C3",   "A1",   "A2",   "A3",   "E1",   "E2",   "E3", "AE1",   "AE2",   "AE3",   "N1",   "N2",   "N3")

#coldata = c("cathepsinG","cathepsinG","hpr3","hpr3","elastase","elastase","N","N")
#replicates = c("C1",    "C2",   "H1",   "H2",    "E1",   "E2", "N1",   "N2")
#treatments = c("AEBSF", "EDTA", "AEBSF+EDTA", "N")

install.packages("BiocManager")
install.packages("edgeR")
install.packages("Logolas")
install.packages("tidyverse")
install.packages("vsn")
install.packages("ggVennDiagram")
install.packages("dplyr")
install.packages("tidyr")
install.packages("expss")
install.packages("EnhancedVolcano")
install.packages("gridExtra")
install.packages("ggdendro")
install.packages("viridis")
install.packages("plotly")
install.packages("factoextra")
install.packages("apeglm")
install.packages("ashr")
install.packages("IHW")
install.packages("PoiClaClu")
install.packages("httr")
install.packages("gtable")
install.packages("munsell")
install.packages("colorspace")
install.packages("scales")
install.packages("vctrs")
install.packages("utf8")
BiocManager::install("DECIPHER")
BiocManager::install("DESeq2")

#load R packages for DEG data analysis
library(DESeq2)
library(apeglm)
library(ashr)
library(IHW)
library(PoiClaClu)
library(edgeR)
library(vctrs)

#load R packages for data visualization
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(plotly)
library(stats)
library(MASS)
library(reshape2)
library(viridis)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(EnhancedVolcano)
library(expss)
library(stringr)
library(tidyverse)
library(vsn)
library(ggVennDiagram)
library(dplyr)
library(Logolas)
library(stringr)

#remove stops and output stop count
data_NS <- subset(data, grepl("X", row.names(data))==FALSE)
orig_ct <- dim(data)[1]
new_ct <- dim(data_NS)[1]
stop_ct <- orig_ct - new_ct
print(paste("Removed",stop_ct, "stop codon substitutions.")); print(paste(round(stop_ct/length(data[,1])*100,digits = 1), "percent of substitutions were stops."))
countData <- as.matrix(data_NS)
colData <- data.frame(condition=factor(coldata))
dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds$condition <- relevel(dds$condition, ref = "N")  #set input counts to be reference level
keep <- rowSums(counts(dds)) >= 10 #rid dds of very-low-count peptides
dds <- dds[keep,]


dds=DESeq(dds)
saveRDS(dds, file = "dds_object.rds")

#The above part could be skipped if the code was ran before and the RDS file was already generated
dds
#This was generated with the above code, but loaded now for speed
dds <- readRDS("C:/Users/enoch/Downloads/Github ALL/dds_object.rds")
dir.create("Significant_Peptides")


for (i in 1:length(treatments)){
  treatments[i]
  res_bonf <- results(dds, contrast=c("condition",treatments[i],"N"), alpha=0.05,pAdjustMethod = "bonferroni")
  pdf(paste0(treatments[i], '_DEseq_plotMA.pdf'))
  plotMA(res_bonf,alpha=.05,ylim=c(-8,8),main='Plot')
  dev.off()
}


for (i in 1:length(treatments)){
  
  print(treatments[i])
  
  #standard FDR 0.05
  res_fdr <- results(dds, contrast=c("condition",treatments[i],"N"), alpha=0.05,pAdjustMethod = "fdr")
  summary(res_fdr)
  res_fdr_sig <- subset(as.data.frame(res_fdr[order(res_fdr$pvalue),]), padj < 0.05)
  AA_seq <- rownames(res_fdr_sig)
  res_fdr_sig  <- cbind(AA_seq, res_fdr_sig)
  write.table(res_fdr_sig,file=paste0("Significant_Peptides/",treatments[i],"_fdr_deseq2.txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  
  #bonferroni 0.05
  res_bonf <- results(dds, contrast=c("condition",treatments[i],"N"), alpha=0.05,pAdjustMethod = "bonferroni")
  summary(res_bonf)
  res_bonf_sig <- subset(as.data.frame(res_bonf[order(res_bonf$pvalue),]), padj < 0.05)
  
  AA_seq <- rownames(res_bonf_sig)
  res_bonf_sig  <- cbind(AA_seq, res_bonf_sig)
  res_bonf_sig_enrich <- subset(as.data.frame(res_bonf_sig[order(res_bonf$baseMean),]), log2FoldChange > 0)
  res_bonf_sig_deplete <- subset(as.data.frame(res_bonf_sig[order(res_bonf$baseMean),]), log2FoldChange < 0)
  write.table(res_bonf_sig,file=paste0("Significant_Peptides/",treatments[i],"_bonf_deseq2.txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  write.table(res_bonf_sig_enrich,file=paste0("Significant_Peptides/","enrich_",treatments[i],"_bonf_deseq2.txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  write.table(res_bonf_sig_deplete,file=paste0("Significant_Peptides/","deplete_",treatments[i],"_bonf_deseq2.txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  
  
  
  #plot most significant gene for each treatment:  
  ix = which.min(res_fdr$padj) # most significant
  png(paste(treatments[i],"_Significant.png", sep=""))
  plotCounts(dds, gene=rownames(dds[ix,]), las=2, intgroup=c("condition"), main=paste("Most significant for",treatments[i],"is",rownames(dds[ix,])))
  dev.off()
  
  #make Venn Diagrams of overlap between p-correction methods
  fdr_sig_AA <- as.character(rownames(res_fdr_sig))
  bonf_sig_AA <- as.character(rownames(res_bonf_sig))

  venn_list <-  list(FDR = fdr_sig_AA,Bonferroni=bonf_sig_AA)
  # plot  
  png(paste(treatments[i],"_FDRvsBonferroni.png", sep=""))
  plot <- ggVennDiagram(venn_list)
  print(plot)
  dev.off()
  
  

}

