#! /usr/bin/Rscript
#YuEnoch 13-11-2022
#PCAprocessor.R 
#Based off code from Dr. Matt Holding

#Purpose: clusters Enriched Counts, as identified from DESeq2, through Principal Component Analysis 

#         Other Plots:
#         Plots the PCA visualization
#         Plots the Goodness of Fit when testing different number of cluster
#         Plots the Loadings of each property at each residue
#Changes across Experiments: alter accordingly
# 1. Mutagenesis for 5 Amino Acids
# 2. NNK Mutagenesis

#Edit the Following: working directory, treatments, file locations
folder  <- 'C:/Users/enoch/Downloads/Github ALL' #<-----------

source(paste0(folder, '/PCA Analysis/SeqSpace main function.R'))
source(paste0(folder, '/PCA Analysis/Plots.R'))
res.prop1 <- readRDS(paste0(folder,'/PCA Analysis/config.rds')) #contains the properties of each Amino Acid
treatments <- c("cathepsinG", "hpr3", "elastase")

#########
# Setup #
#########

res.prop1

colours<-palette(c("blue",       #1
                   "red",        #2
                   "yellow",     #3
                   "white",      #4
                   "darkgrey",   #5
                   "maroon",     #6
                   "orange",     #7
                   "black",      #8
                   "purple",     #9
                   "darkgreen", #10
                   "cyan", #11
                   "bisque", #12
                   "azure", #13
                   "aquamarine", #14
                   "khaki", #15
                   "pink", #16
                   "brown", #17
                   "cadetblue", #18
                   "darkolivegreen", #19
                   "darkorchid4" #20
                   
)) 


##############
# File names #
##############
library(animation)
library(mclust)
library(seqinr)
library(tidyr)
library(rgl)
library(rglwidget)
library(ape)
library(ggplot2)
library(igraph)
library(vctrs)


############
# Analysis #
############

for (i in 1:length(treatments)){

  setwd(folder)  
  name = treatments[i]
  
  SAPCA <- PCA_MSA (MSA = paste('enrich_',name, '_bonf_deseq2.fasta', sep = ''), res.prop = res.prop1, clusterPCs = 3, clusters = 1:10)
  #Conducts the PCA Analysis to cluster the significantly enriched peptides
    #res.prop: contains properties of each Amino Acid (Charge, Disorder, HPATH, RMW)
    #clusterPCs: how many Principal Components to be used
    #clusters: the minimum:maximum number of clusters
  
  newdir <- (paste(folder,"/PCAresults", name, "/",sep="."))
  dir.create(newdir)
  setwd(newdir)
  
  plot_3Dclusters(SAPCA, plotPCs = 1:3, labels = "Consensus") #3D PCA Plot
  rgl.viewpoint( theta = 0, phi = 15, fov = 60, zoom = 0.8)
  rgl.snapshot(paste(name, "_3DPlot_1.png", sep=""), fmt = 'png')
  plot_3Dclusters(SAPCA, plotPCs = 1:3, labels = "Consensus") #3D PCA Plot
  rgl.viewpoint( theta = 90, phi = 90, fov = 60, zoom = 0.8)
  rgl.snapshot(paste(name, "_3DPlot_2.png", sep=""), fmt = 'png')
  plot_3Dclusters(SAPCA, plotPCs = 1:3, labels = "Consensus") #3D PCA Plot
  rgl.viewpoint( theta = 45, phi = 45, fov = 60, zoom = 0.8)
  rgl.snapshot(paste(name, "_3DPlot_3.png", sep=""), fmt = 'png')
  print(name)

  png(paste(name, "_Loadings.png", sep=""))
  plot <- plot_loadings_matrix(SAPCA,PC=1,magnitude=FALSE)
  print(plot) #Loadings of each property at each residue
  dev.off()
  
  png(paste(name, "_GoodnessOfFit.png", sep=""))
  plot_modelfit(SAPCA) #goodness of fit plot
  dev.off()

  write.csv (SAPCA$seq.space.PCA$loadings,           file = paste(name,"Loadings.All.csv",sep="."))
  write.csv (SAPCA$seq.space.PCA$coordinates,        file = paste(name,"Coordinates.csv",sep=".")) 
  clust.sum <- data.frame(SAPCA$seq.space.clusters$classification)
  colnames(clust.sum) <- "cluster" #find way to reduce number of clusters
  clust.sum[["certainty"]] <- apply(SAPCA$seq.space.clusters$likelihoods,1,max)
  write.csv (clust.sum, file = paste(name,"Clusters.csv",sep="."))
  clust.sum$sequences <- row.names(clust.sum)
  
  #Writes the peptides in each cluster into CSV files
  for (i in 1:max(clust.sum[, 1])){
    a = length(clust.sum[, 1])
    csv_fname = paste(name, "_cluster_", i, ".csv", sep = "")
    row <- data.frame("Sequence", "Cluster", "Certainty")
    write.csv(row, file = csv_fname, sep = ",",
              append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
    for (j in 1:a){
      if(clust.sum[j, 1] == i){
        row  <- data.frame(toString(clust.sum[j, 3]), toString(clust.sum[j, 1]), toString(clust.sum[j, 2]))
        write.table(row, file = csv_fname, sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      }
    }
  } 
}
