Analysis of Phage Display Protease and Supernatant Screens
================
Enoch Yu
2024-07-30

# Introduction

A substrate phage display library of randomized 5-mers was constructed,
with each 5-mer having a FLAG peptide at the end. Phages were
immobilized onto anti-FLAG beads, and the phage library was screened
against 3 neutrophil serine proteases: cathepsin G, elastase, and human
proteinase 3 (hPR3). This was followed by 2 controls: with no protease
(baseline phage release) and with a FLAG peptide elution (unselected
phage library). All were conducted in triplicate.

Screening was repeated for the supernatant releasate of activated
neutrophils, in the presence of no inhibitors, AEBSF, EDTA, and
AEBSF+EDTA. The broad-spectrum inhibitor AEBSF inhibits serine
proteases, while EDTA inhibits metalloproteases. The supernatant formed
a complex biological mixture for subsequent analysis.

After incubation, cleaved/released phages are collected, with its DNA
isolated and purified, followed by 2 rounds of PCR and high-throughput
sequencing on Nextseq2000. Sequencing reads were assembled, quality
filtered, and translated into a counts table of cleaved peptides.

<div class="figure" style="text-align: center">

<img src="phagelibraryprep.png" alt="Phage library generation and experimental design" width="40%" height="20%" />
<p class="caption">
Phage library generation and experimental design
</p>

</div>

# Setup for Analysis

This R markdown file will conduct the data analysis for the Deep
Protease Profiling manuscript, using 2 files as the starting point.
These files are the counts of cleaved peptides across proteases and
across supernatant mixtures, generated from the GitHub pipeline with a
row-sum filter of 6. These files were provided as supplementary data
files and can be obtained at our GitHub repository.

neutrophil_mixture_data.csv neutrophil_serine_protease_data.csv

To proceed, please place these 2 files in the same directory as this R
markdown file. This R markdown also should be run after a prior R
markdown file (CombinedAnalysis.Rmd). 6 control frequency files
previously generated would be used in this R markdown file.

### Opening Files

Cleaved peptide counts were passed through a rowsum filter of 6. The
unfiltered counts table was presented in the supplementary data files.

``` r
set.seed(123)
protease_counts <- read.csv(file = "neutrophil_serine_protease_data.csv", header=TRUE,row.names=1)
mixture_counts <- read.csv(file = "neutrophil_mixture_data.csv", header=TRUE,row.names=1)
```

Data Parameters for the protease conditions were loaded: cathepsin G,
elastase, hPR3, no protease (negative control), and FLAG peptide elution
of phage (positive control). Each treatment was conducted in triplicate.

``` r
proteases = c("cathepsin G","elastase","hPR3","no protease")
protease_replicates = c("C1", "C2", "C3", "E1", "E2", "E3", "H1", "H2", "H3", "V1", "V2", "V3", "N1", "N2", "N3")
protease_col = c("cathepsin G","cathepsin G","cathepsin G" ,"elastase","elastase","elastase","hPR3","hPR3","hPR3","no protease", "no protease", "no protease", "FLAG","FLAG","FLAG")
protease_colname = c("cathepsin G 1","cathepsin G 2","cathepsin G 3","elastase 1","elastase 2","elastase 3","hPR3 1","hPR3 2","hPR3 3","no protease 1", "no protease 2", "no protease 3", "FLAG 1","FLAG 2","FLAG 3")
control = "FLAG"
```

Data Parameters for the mixture conditions were loaded: supernatant
mixtures of activated neutrophils with no inhibitor, with AEBSF, with
EDTA, with AEBSF + EDTA, and FLAG peptide elution of phage (positive
control). Each treatment was conducted in triplicate.

``` r
mixtures = c("no inhibitor","AEBSF","EDTA","AEBSF+EDTA")
mixture_replicates = c("NI1","NI2","NI3","A1","A2","A3","ED1","ED2","ED3","AE1","AE2","AE3","G1","G2","G3")
mixture_col = c("no inhibitor","no inhibitor","no inhibitor", "AEBSF","AEBSF","AEBSF","EDTA","EDTA","EDTA","AEBSF+EDTA","AEBSF+EDTA","AEBSF+EDTA","FLAG","FLAG","FLAG")
mixture_colname = c("no inhibitor 1","no inhibitor 2","no inhibitor 3","AEBSF 1","AEBSF 2","AEBSF 3","EDTA 1","EDTA 2","EDTA 3","AEBSF+EDTA 1", "AEBSF+EDTA 2", "AEBSF+EDTA 3", "FLAG 1","FLAG 2","FLAG 3")
```

Necessary packages were loaded: dplyr, tidyr, RColorBrewer, reshape,
factoextra, stringr, openxlsx, stats, multcomp, Hmisc, gridExtra,
tidyverse, ggplot2, ggseqlogo, ggpubr, ggbeeswarm, ggVennDiagram,
ggVenn, pheatmap, PoiClaClu, BiocManager, S4Arrays, DelayedArray,
DESeq2, DECIPHER, plotly, plot3D, dbscan, forcats, umap, immunedeconv,
plotly, clusterSim, immunedeconv

# Initial Analysis of Protease and Mixture Profiles

The cleaved peptides data was analyzed with a bioinformatics approach,
adapting computational packages developed for DNA and RNA sequencing
datasets.

## DESeq2 Analysis Set Up

DESeq2 was employed for data normalization and differential enrichment
analysis to identify peptides significantly cleaved relative to control
(the unselected phage library). Prior to analysis, peptides with a stop
codon or row sum below 10 were removed.

DESeq2 Analysis on Peptide Dataset for Proteases

``` r
data_NS <- subset(protease_counts, grepl("X", row.names(protease_counts))==FALSE)
countData <- as.matrix(data_NS)
colData <- data.frame(condition=factor(protease_col))
dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

``` r
dds$condition <- relevel(dds$condition, ref = "FLAG")
keep <- rowSums(counts(dds)) >= 10
dds_protease <- dds[keep,]
```

DESeq2 Analysis on Peptide Dataset for Mixtures

``` r
data_NS <- subset(mixture_counts, grepl("X", row.names(mixture_counts))==FALSE)
countData <- as.matrix(data_NS)
colData <- data.frame(condition=factor(mixture_col))
dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

``` r
dds$condition <- relevel(dds$condition, ref = "FLAG")
keep <- rowSums(counts(dds)) >= 10
dds_mixture <- dds[keep,]
```

DESeq2 Analysis and File Save for Future Testing

``` r
dds_protease=DESeq(dds_protease)
dds_mixture=DESeq(dds_mixture)
saveRDS(dds_protease, file = paste0("protease_dds_object.rds"))
saveRDS(dds_mixture, file = paste0("mixture_dds_object.rds"))
```

## Clustering Analysis

Poisson distance hierarchal clustering and PCA (principal component
analysis) visualization was conducted on the DESeq2-normalized datasets
to assess intra- and inter-condition variability. Clustering also
provided an overview of broad-scale substrate cleavage patterns,
displaying similarities/differences across conditions.

Poisson Distance Dendrogram for Proteases

``` r
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

poisd_protease <- PoissonDistance(t(counts(dds_protease)))
samplePoisDistMatrix <- as.matrix(poisd_protease$dd)
rownames(samplePoisDistMatrix) <- paste(dds_protease$condition)
colnames(samplePoisDistMatrix) <- NULL
rownames(samplePoisDistMatrix) <- paste(protease_colname)
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd_protease$dd, clustering_distance_cols = poisd_protease$dd, col = colors)
```

Poisson Distance Dendrogram for Mixtures

``` r
poisd_mixture <- PoissonDistance(t(counts(dds_mixture)))
samplePoisDistMatrix <- as.matrix(poisd_mixture$dd)
rownames(samplePoisDistMatrix) <- paste(dds_mixture$condition)
colnames(samplePoisDistMatrix) <- NULL
rownames(samplePoisDistMatrix) <- paste(mixture_colname)
pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd_mixture$dd, clustering_distance_cols = poisd_mixture$dd, col = colors)
```

<div class="figure" style="text-align: center">

<img src="Poisson_protease.png" alt="Poisson Distance Plots for Proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="Poisson_Mixture.png" alt="Poisson Distance Plots for Proteases (Left) and Mixtures (Right)" width="40%" height="20%" />
<p class="caption">
Poisson Distance Plots for Proteases (Left) and Mixtures (Right)
</p>

</div>

PCA Plots for Protease and Mixture conditions

``` r
rld <- rlog(dds_protease, blind=FALSE)
pca <- prcomp(t(assay(rld)))
fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "PCA - Protease Conditions")

rld <- rlog(dds_mixture, blind=FALSE)
pca <- prcomp(t(assay(rld)))
fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "PCA - Mixture Conditions")
```

<div class="figure" style="text-align: center">

<img src="PCA_protease.png" alt="PCA Plots for Proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="PCA_mixture.png" alt="PCA Plots for Proteases (Left) and Mixtures (Right)" width="40%" height="20%" />
<p class="caption">
PCA Plots for Proteases (Left) and Mixtures (Right)
</p>

</div>

# Generation of Significantly Cleaved Peptides

Phage libraries are known to possess baseline compositional bias in
amino acid expression, impacting peptides and sequence counts. To
account for background bias, DESeq2 identified peptides significantly
enriched against the unselected library as control (FLAG peptide
elution). These cleaved peptides were used to capture a proteolytic
activity profile for that protease or mixture.

To limit false positives, a Bonferroni correction with stringent
thresholds of an adjusted p-value \< 10^-9 and log 2-fold change \> 1
were applied. Peptides that met these conditions were classified as a
significantly cleaved peptide.

Defining Parameters and Thresholds

``` r
pvalthreshold = 0.000000001
log2foldthreshold = 1
combined <- append(proteases, mixtures)
peptideList <- list()
```

Differential Enrichment Analysis for Protease Conditions

``` r
for(i in 1:length(proteases)){
  res_bonf <- results(dds_protease, contrast=c("condition",proteases[i],"FLAG"), alpha=pvalthreshold,pAdjustMethod = "bonf")
  res_bonf_sig <- subset(as.data.frame(res_bonf[order(res_bonf$pvalue),]), padj < pvalthreshold)
  AA_seq <- rownames(res_bonf_sig)
  res_bonf_sig  <- cbind(AA_seq, res_bonf_sig)
  res_bonf_sig_enrich <- subset(as.data.frame(res_bonf_sig[order(res_bonf$baseMean),]), log2FoldChange > log2foldthreshold)
  bonf_sig_AA <- as.character(rownames(res_bonf_sig_enrich))
  print(paste0(proteases[i], ": ", length(res_bonf_sig_enrich$baseMean)))
  write.table(res_bonf_sig_enrich,file=paste0("deseq2_", proteases[i],".txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  write.table(bonf_sig_AA,file=paste0("sig_peptides_", proteases[i],".txt"), col.names=FALSE, append=FALSE, quote=FALSE, row.names = FALSE)
  peptideList <- append(peptideList, list(bonf_sig_AA))
  names(peptideList)[length(peptideList)] <- proteases[i]
}
```

    ## [1] "cathepsin G: 5985"
    ## [1] "elastase: 5325"
    ## [1] "hPR3: 1285"
    ## [1] "no protease: 1291"

Differential Enrichment Analysis for Mixture Conditions

``` r
for(i in 1:length(mixtures)){
  res_bonf <- results(dds_mixture, contrast=c("condition",mixtures[i],"FLAG"), alpha=pvalthreshold,pAdjustMethod = "bonf")
  res_bonf_sig <- subset(as.data.frame(res_bonf[order(res_bonf$pvalue),]), padj < pvalthreshold)
  AA_seq <- rownames(res_bonf_sig)
  res_bonf_sig  <- cbind(AA_seq, res_bonf_sig)
  res_bonf_sig_enrich <- subset(as.data.frame(res_bonf_sig[order(res_bonf$baseMean),]), log2FoldChange > log2foldthreshold)
  bonf_sig_AA <- as.character(rownames(res_bonf_sig_enrich))
  print(paste0(mixtures[i], ": ", length(res_bonf_sig_enrich$baseMean)))
  write.table(res_bonf_sig_enrich,file=paste0("deseq2_", mixtures[i],".txt"), col.names=TRUE, append=FALSE, quote=FALSE, row.names = FALSE)
  write.table(bonf_sig_AA,file=paste0("sig_peptides_", mixtures[i],".txt"), col.names=FALSE, append=FALSE, quote=FALSE, row.names = FALSE)
  peptideList <- append(peptideList, list(bonf_sig_AA))
  names(peptideList)[length(peptideList)] <- mixtures[i]
}
```

    ## [1] "no inhibitor: 5166"
    ## [1] "AEBSF: 558"
    ## [1] "EDTA: 5190"
    ## [1] "AEBSF+EDTA: 313"

MA plots

To visualize the peptides and their corresponding counts or log 2 fold
change, relative to control, MA plots were made. Blue dots represent
peptides with a Bonferroni p-adj \< 10^-9.

``` r
for (i in 1:length(proteases)){
  res_bonf <- results(dds_protease, contrast=c("condition",proteases[i],control), alpha=pvalthreshold,pAdjustMethod = "bonferroni")
  pdf(paste0(proteases[i], '_DEseq_plotMA.pdf'))
  plotMA(res_bonf,alpha=pvalthreshold,ylim=c(-8,8),main='Plot')
  dev.off()
}
for (i in 1:length(mixtures)){
  res_bonf <- results(dds_mixture, contrast=c("condition",mixtures[i],control), alpha=pvalthreshold,pAdjustMethod = "bonferroni")
  pdf(paste0(mixtures[i], '_DEseq_plotMA.pdf'))
  plotMA(res_bonf,alpha=pvalthreshold,ylim=c(-8,8),main='Plot')
  dev.off()
}
```

<div class="figure" style="text-align: center">

<img src="csgmaplot.png" alt="Example of a MA plot for cathepsin G" width="40%" height="20%" />
<p class="caption">
Example of a MA plot for cathepsin G
</p>

</div>

Sequence Logos

For the significantly cleaved peptides, sequence logos were constructed
using ggseqlogo to depict amino acid motifs.

``` r
proteases_list <- list(`cathepsin G` = peptideList$`cathepsin G`,elastase=peptideList$elastase, hPR3=peptideList$hPR3, `no protease`=peptideList$`no protease`)
mixtures_list <- list(`no inhibitor` = peptideList$`no inhibitor`,AEBSF=peptideList$AEBSF, EDTA=peptideList$EDTA, `AEBSF+EDTA`=peptideList$`AEBSF+EDTA`)
combined_list <- append(proteases_list, mixtures_list)
ggseqlogo(combined_list)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-10-1.png" width="60%" />

``` r
ggseqlogo(combined_list, method = 'prob')
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-10-2.png" width="60%" />

## Peptide Space Comparisons

Venn diagrams using ggVennDiagram were used to determine the overlap in
cleaved peptides between proteases, mixtures, and other conditions, as
well as areas of distinct activity.

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-11-1.png" width="49%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-11-2.png" width="49%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-11-3.png" width="49%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-11-4.png" width="49%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-11-5.png" width="49%" />

Functions to Convert Peptides between TXT and Fasta files

``` r
txt_to_fasta <- function(fname) {
  AA_seqs <- scan(file=paste0(fname, ".txt"), what="", sep="\n")
  fastaFile <- paste0(fname,".fasta")
  fileConn<-file(fastaFile)
  sink(fastaFile)
  AA_seqs <- as.data.frame(AA_seqs)
  for (i in 1:length(AA_seqs$AA_seqs)){
    cat(paste0(">",AA_seqs$AA_seqs[i]))
    cat("\n")
    cat(AA_seqs$AA_seqs[i])
    cat("\n")
  }
  close(fileConn)
  sink() 
}

fasta_to_txt <- function(fname) {
  AA_seqs <- scan(file=paste0(fname, ".fasta"), what="", sep="\n")
  AA_seqs <- AA_seqs[!grepl(">", AA_seqs)]
  txtFile <- paste0(fname,".txt")
  fileConn<-file(txtFile)
  sink(txtFile)
  AA_seqs <- as.data.frame(AA_seqs)
  for (i in 1:length(AA_seqs$AA_seqs)){
    cat(AA_seqs$AA_seqs[i])
    cat("\n")
  }
  close(fileConn)
  sink() 
}
```

## Unsupervised Alignment of Significantly Cleaved Peptides, using Decipher

In this 5-mer peptide format, a key limitation was that proteolytic
cleavage may have occurred anywhere along the 5-mer. Based on 4 possible
P1-P1’ cleavage sites, this formed a 9-mer alignment space. The DECIPHER
R Package was used to perform unsupervised alignment and estimate the P1
position.

This alignment was conducted on the significantly cleaved peptides
previously determined for each condition. The parameter terminalGap was
adjusted to fit peptides into 9-mer space.

``` r
for(i in 1:length(combined)){
  txt_to_fasta(paste0("sig_peptides_",combined[i]))
  fastaFile <- paste0("sig_peptides_",combined[i], ".fasta")
  seqs <- readAAStringSet(fastaFile)
  aligned <- AlignSeqs(seqs, gapOpening = -200, useStructures = TRUE, terminalGap = -8)
  #the parameter terminalGap was adjusted to fit 9-mer space
  writeXStringSet(aligned,file= paste0("aligned_sig_peptides_",combined[i], ".fasta"))
  fasta_to_txt(paste0("aligned_sig_peptides_",combined[i]))
}
```

Sequence Logos of Aligned Peptides

The aligned position with the most dominant motifs was determined as P1
for the 9-mer sequence.

``` r
alignedpeptideList <- list()
for(i in 1:length(combined)){
  alignedAA <- read.table(paste0("aligned_sig_peptides_",combined[i], ".txt"))
  df_aligned <- apply(alignedAA, 1, function(row) {
    unlist(strsplit(row, ""))
  })
  df_aligned <- as.data.frame(t(df_aligned))

  seqs <- readAAStringSet(paste0("aligned_sig_peptides_",combined[i], ".fasta"))
  consensus <- consensusMatrix(seqs, as.prob = T)
  consensus <- consensus[rownames(consensus) != "-", ]
  highest_col <- which.max(colSums(consensus * (consensus > 0.095)))
  df_9mer <- df_aligned[, (highest_col-4):(highest_col+4)]
  merged_9mer <- apply(df_9mer, 1, function(row) paste(row, collapse = ""))
  alignedpeptideList <- append(alignedpeptideList, list(merged_9mer))
  names(alignedpeptideList)[length(alignedpeptideList)] <- combined[i]
  write.table(merged_9mer, file = paste0("aligned_9mer_",combined[i], ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  txt_to_fasta(paste0("aligned_9mer_",combined[i]))
}
```

Logos of aligned peptides were graphed using ggseqlogo.

``` r
proteases_aligned_list <- list(`cathepsin G` = alignedpeptideList$`cathepsin G`,elastase=alignedpeptideList$elastase, hPR3=alignedpeptideList$hPR3, `no protease`=alignedpeptideList$`no protease`)
mixtures_aligned_list <- list(`no inhibitor` = alignedpeptideList$`no inhibitor`,AEBSF=alignedpeptideList$AEBSF, EDTA=alignedpeptideList$EDTA, `AEBSF+EDTA`=alignedpeptideList$`AEBSF+EDTA`)
combined_aligned_list <- append(proteases_aligned_list, mixtures_aligned_list)
ggseqlogo(combined_aligned_list)
```

![](Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggseqlogo(combined_aligned_list, method = 'prob')
```

![](Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

## Synthetic Alignment of Significantly Cleaved Peptides

An alternative supervised alignment method was synthetic alignment,
centering peptides at the P1 position based on each protease’s known P1
site preferences. This was conducted for the 3 neutrophil serine
proteases. Peptides without the specified amino acids were unaligned.

Cathepsin G: F/Y Elastase: V/I hPR3: V/I

``` r
strfind <- function(s, char) { #function for finding the first index where an amino acid appears in the peptide
  temp <- 1
  if(length(rev(grep(char, unlist(strsplit(s, NULL)), fixed=T))) > 1){
    temp <- length(rev(grep(char, unlist(strsplit(s, NULL)), fixed=T)))}
  rev(grep(char, unlist(strsplit(s, NULL)), fixed=T))[temp]}
```

Synthetic Alignment

``` r
centers <- list(c("F", "Y"), c("V", "I"), c("V", "I"))
syntheticList <- list()
unalignedList <- list()
for(i in 1:3){
  sequences <- peptideList[[i]]
  synthetic <- c()
  unaligned <- c()
  for(j in 1:length(sequences)){
    pep = sequences[j]
    for(k in 1:length(centers)){
      index <- strfind(pep, centers[[i]][k])
      if(!is.na(index)){break}}
    if(is.na(index)){unaligned <- c(unaligned, pep)
    next}
    left <- 9 - nchar(pep) - index + 1
    right <- index - 1
    alignedpep <- paste0(paste0(rep('-', left), collapse = ""), pep, paste0(rep('-', right), collapse = ""))
    synthetic <- c(synthetic, alignedpep)}
  syntheticList <- append(syntheticList, list(synthetic))
  names(syntheticList)[length(syntheticList)] <- combined[i]
  unalignedList <- append(unalignedList, list(unaligned))
  names(unalignedList)[length(unalignedList)] <- combined[i]
  write.table(synthetic, file = paste0("synthetic_9mer_",combined[i], ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  txt_to_fasta(paste0("synthetic_9mer_",combined[i]))
}
```

Sequence logos of Synthetically Aligned Peptides

``` r
synthetic_aligned_list <- list(`cathepsin G` = syntheticList$`cathepsin G`,elastase=syntheticList$elastase, hPR3=syntheticList$hPR3)
synthetic_unaligned_list <- list(`cathepsin G` = unalignedList$`cathepsin G`,elastase=unalignedList$elastase, hPR3=unalignedList$hPR3)

ggseqlogo(synthetic_aligned_list) #aligned peptides
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-18-1.png" width="30%" />

``` r
ggseqlogo(synthetic_unaligned_list) #unaligned peptides
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-18-2.png" width="30%" />

## Icelogos for Peptides with Unsupervised Alignment, referencing the background unselected library

To assess the relative enrichment of peptides at each position, compared
to the background phage library, Icelogos were further generated using
the Icelogo SOAP server. Motifs with a p\<0.05 for % difference,
relative to control

<https://github.com/compomics/icelogo>

Positive Set: significantly cleaved peptides with DECIPHER alignment
Negative Set: list of peptides representing amino acid frequencies of
the unselected phage library (FLAG peptide elution)

Fasta files of the Negative Set was sourced from the Preliminary
Analysis Rmd file.

Control Frequencies for Protease Screen:
‘control_frequency_protease.fasta’ Control Frequencies for Mixture
Screen: ‘control_frequency_mixture.fasta’ Control Frequencies across
both Screens: ‘control_frequency_combined.fasta’

<div class="figure" style="text-align: center">

<img src="cathepsinG_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="elastase_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="hPR3_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="no protease_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="noinhibitor_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="AEBSF_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="EDTA_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="AEBSFEDTA_icelogo-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" />
<p class="caption">
PCA Plots for proteases (Left) and Mixtures (Right)
</p>

</div>

## Icelogos for Peptides with Synthetic Alignment

For 3 proteases with synthetic alignment, to assess if the screen
provided a true signal, a similarly synthetically aligned background
dataset was used. Every possible peptide from the unselected phage
library was synthetically aligned based on F/Y or V/I. This dataset was
used as the Negative Set for Icelogo.

``` r
data_NS <- subset(protease_counts, grepl("X", row.names(protease_counts))==FALSE)
rowsum_filter <- data_NS[rowSums(data_NS) >= 10, ]
FLAG_peptides <- rownames(rowsum_filter[rowSums(rowsum_filter[, 13:15]) > 0, ])
sequences <- FLAG_peptides

centers <- list(c("F", "Y"), c("V", "I"))
syntheticFLAG <- list()
unalignedFLAG <- list()
for(i in 1:2){
  synthetic <- c()
  unaligned <- c()
  for(j in 1:length(sequences)){
    pep = sequences[j]
    for(k in 1:length(centers)){
      index <- strfind(pep, centers[[i]][k])
      if(!is.na(index)){break}}
    if(is.na(index)){unaligned <- c(unaligned, pep)
    next}
    left <- 9 - nchar(pep) - index + 1
    right <- index - 1
    alignedpep <- paste0(paste0(rep('-', left), collapse = ""), pep, paste0(rep('-', right), collapse = ""))
    synthetic <- c(synthetic, alignedpep)}
  syntheticFLAG <- append(syntheticFLAG, list(synthetic))
  names(syntheticFLAG)[length(syntheticFLAG)] <- i
  unalignedFLAG <- append(unalignedFLAG, list(unaligned))
  names(unalignedFLAG)[length(unalignedFLAG)] <- i
  write.table(synthetic, file = paste0("synthetic_9mer_FLAG", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  txt_to_fasta(paste0("synthetic_9mer_FLAG", i,combined[i]))
  write.table(unaligned, file = paste0("unaligned_synthetic_9mer_FLAG", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  txt_to_fasta(paste0("unaligned_synthetic_9mer_FLAG", i,combined[i]))
}

write.table(synthetic, file = paste0("synthetic_9mer_FLAG", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(unaligned, file = paste0("unaligned_synthetic_9mer_FLAG", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
txt_to_fasta(paste0("synthetic_9mer_FLAG", i))
txt_to_fasta(paste0("unaligned_synthetic_9mer_FLAG", i))
```

<div class="figure" style="text-align: center">

<img src="cathepsinG_synthetic-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="elastase_synthetic-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" /><img src="hPR3_synthetic-1.png" alt="PCA Plots for proteases (Left) and Mixtures (Right)" width="40%" height="20%" />
<p class="caption">
PCA Plots for proteases (Left) and Mixtures (Right)
</p>

</div>

# Motif Analysis of Substrate Peptides from Mixtures

To analyze distinct patterns of proteolytic activity, the significantly
cleaved peptides of supernatant mixtures were assessed with
dimensionalty reduction techniques PCA and UMAP. Peptides were mapped by
each position based on 4 amino acid characteristics (charge, disorder,
hpath, and relative molecular weight).

PCA/UMAP was followed by K-Means clustering, with the average silouette
width used to estimate the optimal number of clusters.

Defining Amino Acid Characteristics

``` r
AAreference <- data.frame(
  row.names = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "-"),
  CHRG = c(0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, NA, NA),
  DISORD = c(0.0600, 0.0200, 0.1920, 0.7360, -0.6970, 0.1660, 0.3030, -0.4860, 0.5860, -0.3260, -0.3970, 0.0070, 0.9870, 0.3180, 0.1800, 0.3410, 0.0590, -0.1210, -0.8840, -0.5100, NA, NA),
  HPATH = c(1.80, 2.50, -3.50, -3.50, 2.80, -0.40, -3.20, 4.50, -3.90, 3.80, 1.90, -3.50,  -1.60, -3.50, -4.50, -0.80, -0.70, 4.20, -0.90, -1.30, NA, NA),
  RMW = c(15.090, 47.160, 59.100, 73.130, 91.190, 1.070, 81.160, 57.180, 72.190, 57.180, 75.210, 58.120, 41.130, 72.150, 100.200, 31.090, 45.120, 43.150, 130.230, 107.190, NA, NA),
  NOTGAP = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0))
properties_df <- as.data.frame(AAreference)[, -ncol(AAreference)]
```

The subsequent code, from dataframe generation to ggseqlogo creation,
was repeated for each of the 4 mixture conditions with different
inhibitors.

## Generation of Input Dataframe

A 20 column dataframe, of 4 amino acid properties across 5 positions for
each peptide, was generated. This was normalized based on columns and
used as input for PCA and UMAP.

``` r
i = 5 #i = 5-8 for mixtures
peptides = peptideList[i]

features_df <- data.frame(matrix(nrow = length(peptides), ncol = nchar(peptides[[1]][1])*ncol(properties_df)))
colnames(features_df) <- rep(paste0(colnames(properties_df), rep(1:nchar(peptides[[1]][1]), each = ncol(properties_df))), length(peptides))
for (k in 1:length(peptides[[1]])) {
  rows <- data.frame(col1 = NA)
  for (j in 1:nchar(peptides[[1]][1])){
    amino_acid <- substr(peptides[[1]][k], j, j)
    rows_with <- properties_df[rownames(properties_df) == amino_acid, ]
    rows <- cbind(rows, rows_with)
  }
  rows <- rows[,-1]
  colnames(rows) <- rep(paste0(colnames(properties_df), rep(1:nchar(peptides[[1]][1]), each = ncol(properties_df))), length(peptides))
  rownames(rows) <- c(peptides[[1]][k])
  features_df <- rbind(features_df, rows)
}
features_df <- features_df[-1, ]
df.n <- clusterSim::data.Normalization(features_df, type="n1", normalization="column")
```

## PCA Analysis

After silhouette Width analysis, the number of clusters for K-means
clustering was adjusted. An example plot of cleaved peptides by the
supernatant without inhibitors was presented, with all plots in the
supplemental figures.

``` r
pca_result <- prcomp(df.n, scale. = TRUE)
fviz_nbclust(pca_result$x, FUNcluster=kmeans, k.max = 20, method = "silhouette")
kmeans_result <- kmeans(pca_result$x, centers = 9) #number of centrers adjusted based on the peak average silhouette width
df.n$PCACluster <- kmeans_result$cluster
pca_result$cluster <- as.factor(kmeans_result$cluster)

plot_pca <- plot_ly(x = pca_result$x[,1], y = pca_result$x[,2], type = 'scatter', mode = 'markers',color = pca_result$cluster, colors = rainbow(5), marker = list(size = 5))
plot_pca
```

<div class="figure" style="text-align: center">

<img src="mixturepca1.png" alt="PCA Visualization of significantly cleaved peptides by the supernatant of activated neutrophils" width="30%" height="20%" /><img src="mixturepca2.png" alt="PCA Visualization of significantly cleaved peptides by the supernatant of activated neutrophils" width="30%" height="20%" />
<p class="caption">
PCA Visualization of significantly cleaved peptides by the supernatant
of activated neutrophils
</p>

</div>

## UMAP Analysis

The prior analysis was repeated for UMAP. All plots were placed in the
supplemental figures.

``` r
umap_result <- umap(df.n, scale = TRUE)
fviz_nbclust(umap_result$layout, FUNcluster=kmeans, k.max = 20, method = "silhouette")
kmeans_result <- kmeans(umap_result$layout, centers = 4) #number of centrers adjusted based on the peak average silhouette width
df.n$UMAPCluster <- kmeans_result$cluster
umap_result$cluster <- as.factor(kmeans_result$cluster)

plot_umap <- plot_ly(x = umap_result$layout[,1], y = umap_result$layout[,2], type = 'scatter', mode = 'markers',color = umap_result$cluster, colors = rainbow(5), marker = list(size = 5))
plot_umap
```

<div class="figure" style="text-align: center">

<img src="mixtureumap1.png" alt="UMAP Visualization of Significantly Cleaved Peptides by the supernatant of activated neutrophils" width="30%" height="20%" /><img src="mixtureumap2.png" alt="UMAP Visualization of Significantly Cleaved Peptides by the supernatant of activated neutrophils" width="30%" height="20%" />
<p class="caption">
UMAP Visualization of Significantly Cleaved Peptides by the supernatant
of activated neutrophils
</p>

</div>

## Assessing Motif Logos

The sequence logos for each cluster from PCA and UMAP analysis were
generated

``` r
temp <- df.n
temp$Peptide <- row.names(df.n)
PCAclusters <- split(temp$Peptide, temp$PCACluster)
names(PCAclusters) <- sapply(1:length(PCAclusters), function(i) {
  paste0(names(PCAclusters)[i], " (", length(PCAclusters[[i]]), ")")
})

PCApeptides <- unlist(PCAclusters)
clustercol <- sub("\\s.*", "", rep(names(PCAclusters), lengths(PCAclusters)))
PCAdf <- data.frame(Peptide = PCApeptides, Cluster = clustercol)
write.table(PCAdf, file = paste0(combined[i],"_PCAclusters.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
```

``` r
temp <- df.n
temp$Peptide <- row.names(df.n)
UMAPclusters <- split(temp$Peptide, temp$UMAPCluster)
names(UMAPclusters) <- sapply(1:length(UMAPclusters), function(i) {
  paste0(names(UMAPclusters)[i], " (", length(UMAPclusters[[i]]), ")")
})

UMAPpeptides <- unlist(UMAPclusters)
clustercol <- sub("\\s.*", "", rep(names(UMAPclusters), lengths(UMAPclusters)))
UMAPdf <- data.frame(Peptide = UMAPpeptides, Cluster = clustercol)
write.table(UMAPdf, file = paste0(combined[i],"_UMAPclusters.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
```

``` r
ggseqlogo(PCAclusters, ncol = 5)
ggseqlogo(UMAPclusters, ncol = 5)
```

<div class="figure" style="text-align: center">

<img src="mixture_pcalogos.png" alt="Clustered Motifs of Significantly Cleaved Peptides by the supernatant, from PCA (Left) and UMAP (Right) with K-means clustering" width="30%" height="20%" /><img src="mixture_umaplogos.png" alt="Clustered Motifs of Significantly Cleaved Peptides by the supernatant, from PCA (Left) and UMAP (Right) with K-means clustering" width="30%" height="20%" />
<p class="caption">
Clustered Motifs of Significantly Cleaved Peptides by the supernatant,
from PCA (Left) and UMAP (Right) with K-means clustering
</p>

</div>

## Analysis between Proteases

In an attempt to detect individual protease signatures within the
complex mixture sample, dimensionalty reduction of PCA/UMAP was applied
on the substrates of two sources, 1) pooled cleaved peptides of all
three serine proteases, and 2) the mixture supernatant. Peptides from
each of the three proteases were labelled by color.

### Pooled Peptides of All Serine Proteases

``` r
df <- enframe(peptideList[1:3], name = "protease", value = "peptide") %>%
  unnest_longer(peptide)
features_df <- data.frame(matrix(nrow = 1, ncol = 5*ncol(properties_df)))
colnames(features_df) <- rep(paste0(colnames(properties_df), rep(1:5, each = ncol(properties_df))), 1)

for (k in 1:nrow(df)) {
  rows <- data.frame(col1 = NA)
  for (j in 1:5){
    amino_acid <- substr(df[k,2], j, j)
    rows_with <- properties_df[rownames(properties_df) == amino_acid, ]
    rows <- cbind(rows, rows_with)
  }
  rows <- rows[,-1]
  colnames(rows) <- rep(paste0(colnames(properties_df), rep(1:5, each = ncol(properties_df))), 1)
  rownames(rows) <- c(df[k,2])
  features_df <- rbind(features_df, rows)
}
features_df <- features_df[-1, ]
df.n <- clusterSim::data.Normalization(features_df, type="n1", normalization="column")
df.n <- bind_cols(df.n, df$protease)
colnames(df.n)[ncol(df.n)] <- "protease"

pca_result <- prcomp(df.n[-ncol(df.n)], scale. = TRUE)
umap_result <- umap(df.n[-ncol(df.n)], scale = TRUE)

plot_pca <- plot_ly(x = pca_result$x[,1], y = pca_result$x[,2], type = 'scatter', mode = 'markers',color = df.n$protease, colors = rainbow(3), marker = list(size = 5))
plot_pca
plot_umap <- plot_ly(x = umap_result$layout[,1], y = umap_result$layout[,2], type = 'scatter', mode = 'markers',color = df.n$protease, colors = rainbow(3), marker = list(size = 5))
plot_umap
```

<img src="combinedpcapeptides.png" width="30%" height="20" style="display: block; margin: auto;" /><img src="combinedumappeptides.png" width="30%" height="20" style="display: block; margin: auto;" />

### Serine Proteases within Cleaved Peptides by Neutrophil Supernatant

``` r
cathepsinG_peptides <- intersect(peptides, peptideList[[1]])
elastase_peptides <- intersect(peptides, peptideList[[2]])
hPR3_peptides <- intersect(peptides, peptideList[[3]])
other_peptides <- setdiff(peptides, union(peptideList[[1]], union(peptideList[[2]], peptideList[[3]])))
mixture_list <- list(cathepsinG_peptides, elastase_peptides, hPR3_peptides, other_peptides)
names(mixture_list) <- c("cathepsin G", "elastase", "hPR3", "other peptides")
df <- enframe(mixture_list, name = "protease", value = "peptide") %>%
  unnest_longer(peptide)

features_df <- data.frame(matrix(nrow = 1, ncol = 5*ncol(properties_df)))
colnames(features_df) <- rep(paste0(colnames(properties_df), rep(1:5, each = ncol(properties_df))), 1)

for (k in 1:nrow(df)) {
  rows <- data.frame(col1 = NA)
  for (j in 1:5){
    amino_acid <- substr(df[k,2], j, j)
    rows_with <- properties_df[rownames(properties_df) == amino_acid, ]
    rows <- cbind(rows, rows_with)
  }
  rows <- rows[,-1]
  colnames(rows) <- rep(paste0(colnames(properties_df), rep(1:5, each = ncol(properties_df))), 1)
  rownames(rows) <- c(df[k,2])
  features_df <- rbind(features_df, rows)
}
features_df <- features_df[-1, ]
df.n <- clusterSim::data.Normalization(features_df, type="n1", normalization="column")
df.n <- bind_cols(df.n, df$protease)
colnames(df.n)[ncol(df.n)] <- "protease"

pca_result <- prcomp(df.n[-ncol(df.n)], scale. = TRUE)
umap_result <- umap(df.n[-ncol(df.n)], scale = TRUE)

plot_pca <- plot_ly(x = pca_result$x[,1], y = pca_result$x[,2], type = 'scatter', mode = 'markers',color = df.n$protease, colors = c("red", "green", "blue", 'yellow'), marker = list(size = 5))
plot_pca
plot_umap <- plot_ly(x = umap_result$layout[,1], y = umap_result$layout[,2], type = 'scatter', mode = 'markers',color = df.n$protease, colors = c("red", "green", "blue", 'yellow'), marker = list(size = 5))
plot_umap
```

<img src="mixturepcapeptides.png" width="30%" height="20" style="display: block; margin: auto;" /><img src="mixtureumappeptides.png" width="30%" height="20" style="display: block; margin: auto;" />

# Physiological Substrate Validation and Prediction

Based on the derived proteolytic activity profile from significantly
cleaved peptides, potential physiological substrates were determined in
silico. Alphafold2 structures of ~18000 unique human protein products
were used. The dssp program calculated secondary structure and solvent
accessibility estimates for each amino acid position.

The total accessibility score for each 5-mer in the protein was
multiplied with the positional weight matrix of the 5-mer substrate
profile, giving a ‘cut-score’ for each protein.

To validate this approach, the MEROPs database of known substrates
categorizes the cleavability of each protein by a specific protease.

## Generation of a Positional Weight Matrix (PWM)

The Positional Weight Matrix represented the frequencies of every amino
acid in each position, in log2 space and relative to control (the
unselected phage library). The frequencies of the FLAG control condition
were derived from the Preliminary Analysis R markdown file. This PWM was
used for subsequent validation and prediction of physiological
substrates.

Control Frequencies for Protease Screen:
‘control_frequency_protease.txt’ Control Frequencies for Mixture Screen:
‘control_frequency_mixture.txt’ Control Frequencies across both Screens:
‘control_frequency_combined.txt’

### Obtaining Amino Acid Frequency data for the Control

``` r
AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

control_frequency <- read.table("control_frequency_protease.txt", header = TRUE)
rownames(control_frequency) <- control_frequency[, 1]
control_frequency <- control_frequency[, -1]
control_frequency <- control_frequency[rownames(control_frequency)!="X",]
control_frequency <- control_frequency[, -ncol(control_frequency)]
protease_control_prob <- control_frequency/colSums(control_frequency)

control_frequency <- read.table("control_frequency_mixture.txt", header = TRUE)
rownames(control_frequency) <- control_frequency[, 1]
control_frequency <- control_frequency[, -1]
control_frequency <- control_frequency[rownames(control_frequency)!="X",]
control_frequency <- control_frequency[, -ncol(control_frequency)]
mixture_control_prob <- control_frequency/colSums(control_frequency)
```

### Generation of PWM for Significantly Cleaved Peptides (unaligned, 5-mers)

``` r
for (i in 1:length(combined)){
  sequences = peptideList[[i]]
  if(combined[i] %in% unlist(proteases)){
    control_prob <- protease_control_prob}
  else{control_prob <- mixture_control_prob}
  
  AA_counts <- as.data.frame(matrix(nrow=nchar(length(sequences[1])), ncol=20))
  colnames(AA_counts)<-AA
  for (j in 1:nchar(sequences[i])) {
    temp <- substr(sequences, j, j)
    peptideNames <- AA %>%  unique() %>%data.frame(stringsAsFactors = F)
    colnames(peptideNames) <- "peptideNames"
    result<-apply(peptideNames,1,function(x) str_count(temp, x["peptideNames"]))
    colnames(result) <- peptideNames$peptideNames
    result<-result[, AA]
    for (k in 1:20) {
      AA_counts[j,k] <- sum(result[,k])
    }
  }
  AA_counts$total <- apply(AA_counts, 1, sum)
  PPM <- AA_counts/AA_counts$total #position probability matrix 
  PPM <- PPM[,1:20]
  PPM <- PPM %>% t() %>% as.data.frame()
  colnames(PPM) <- as.character(1:nchar(sequences[1]))

  pwm <- log2((PPM/control_prob))
  colnames(pwm) <- as.character(1:nchar(sequences[1]))
  pwm_output <- t(pwm)
  pwm_output <- replace(pwm_output, is.infinite(pwm_output), 0)
  write.table(pwm_output, file = paste0(combined[i], '_pwm.txt'), sep = "\t", row.names = FALSE)
}
```

Example PWM for AEBSF+EDTA Mixture Condition

``` r
pwm_output
```

    ##            A          C          D          E          F          G          H
    ## 1 -0.1995469 -1.8618035 -0.9722347 -2.2106657  0.7289047 -1.1758030 -0.4955468
    ## 2 -0.0839636  0.3658864 -2.4429397 -0.3255187  0.7215508 -1.2267292 -0.7868103
    ## 3  0.5001332  0.6681490 -3.5353484 -3.2079338 -0.2039981 -2.1622572 -1.7242127
    ## 4  0.3333831  0.1162614 -3.3880542 -1.3239097  0.8098335 -2.3700648 -2.2881587
    ## 5  0.2171522 -0.6412229 -0.8473906 -2.9746375  0.9146838 -0.7062996 -1.1752196
    ##            I          K          L           M          N         P          Q
    ## 1  0.4119210 0.83133886  0.1897784  0.21062365 -0.3216139 -2.501886 -1.1088586
    ## 2  0.2103042 0.99340063  0.5799720 -1.06604304 -0.7973798 -1.911750 -0.9622788
    ## 3 -0.2546420 1.27472134  0.4251034  0.16506781 -1.7547883 -3.510340 -1.8600432
    ## 4  0.2301400 0.80259142 -0.8570150 -0.09061232 -2.6727604 -1.269944 -1.6110125
    ## 5  0.5290625 0.03023376  0.3643589  0.07731044 -0.6398946 -2.805899 -1.6204591
    ##           R          S          T          V          W          Y
    ## 1 0.7469068 -0.8149775 -1.2768802  0.1432961  1.1375137  0.9421939
    ## 2 1.1886381 -0.9527802 -0.8151200 -0.3691384  0.4178152  0.4454659
    ## 3 1.9639667 -1.1959496 -2.0199917  0.6051118 -2.2616464 -0.4536381
    ## 4 1.7013448 -1.2837178 -1.0816984  0.2884817  0.5474301  0.7370825
    ## 5 1.1105180 -0.8651160 -0.6764551  0.2232195  0.7217041  0.2499968

### Generation of PWM for Significantly Cleaved Peptides (aligned with DECIPHER, 9-mers)

``` r
protease_9control_prob <- matrix(rep(rowMeans(protease_control_prob, na.rm = TRUE), each = 9), ncol = 9, byrow = TRUE)
mixture_9control_prob <- matrix(rep(rowMeans(mixture_control_prob, na.rm = TRUE), each = 9), ncol = 9, byrow = TRUE)

for (i in 1:length(combined)){
  sequences = alignedpeptideList[[i]]
  if(combined[i] %in% unlist(proteases)){
    control_prob <- protease_9control_prob
  }
  else{control_prob <- mixture_9control_prob}
  AA_counts <- as.data.frame(matrix(nrow=nchar(length(sequences[1])), ncol=20))
  colnames(AA_counts)<-AA
  for (j in 1:nchar(sequences[i])) {
    temp <- substr(sequences, j, j)
    peptideNames <- AA %>%  unique() %>%data.frame(stringsAsFactors = F)
    colnames(peptideNames) <- "peptideNames"
    result<-apply(peptideNames,1,function(x) str_count(temp, x["peptideNames"]))
    colnames(result) <- peptideNames$peptideNames
    result<-result[, AA]
    for (k in 1:20) {
      AA_counts[j,k] <- sum(result[,k])
    }
  }
  AA_counts$total <- apply(AA_counts, 1, sum)
  PPM <- AA_counts/AA_counts$total
  PPM <- PPM[,1:20]
  PPM <- PPM %>% t() %>% as.data.frame()
  colnames(PPM) <- as.character(1:nchar(sequences[1]))

  pwm <- log2((PPM/control_prob))
  colnames(pwm) <- as.character(1:nchar(sequences[1]))
  pwm_output <- t(pwm)
  pwm_output <- replace(pwm_output, is.infinite(pwm_output), 0)
  write.table(pwm_output, file = paste0(combined[i], '_aligned_pwm.txt'), sep = "\t", row.names = FALSE)
}
```

## Alphafold2 Testing of Proteins

Applying this generated PWM,

``` r
#Alphafold Script etc.
```

# Protease Mixture Deconvolution with CibersortX and EPIC

Based on prior analyses, peptide profiles were identified to broadly
capture the proteolytic signatures of different proteases and mixture
conditions. However, major challenges remain in deconvolution and
differentiating proteases within mixtures objectively.

Drawing upon developed Bulk RNA-seq Deconvolution methods, the peptide
dataset of proteases and mixtures were applied towards this purpose. To
ensure robustness, 2 different algorithms CibersortX and EPIC, were
used.

<https://cibersortx.stanford.edu/> <https://epic.gfellerlab.org/>

A comprehensive protease reference profile was constructed, using the
normalized counts of all cleaved peptides in each protease condition.
This profile was directly applied towards analyzing the separate mixture
dataset of activated neutrophils, determining estimates of proteolytic
activity levels for each protease.

Reference Profile 1) 3 neutrophil serine proteases, 2 phage controls,
ADAMTS13, thrombin

The protease profile consisted of cleaved/released phages in 7
conditions: 3 neutrophil serine proteases (cathepsin G, elastase, hPR3),
no protease (baseline phage release), FLAG peptide elution (unselected
phage library), and 2 blood proteases (ADAMTS13, thrombin).

The 2 phage controls were included as a reference for phage library
background and noise. The 2 blood proteases were included to provide a
control against false signals, since both have no expected
activity/presence in neutrophil supernatants.

Data for ADAMTS13 and thrombin were sourced from a 6-mer NNK phage
display library with a similar composition by Kretz et al., 2018, with
6-mers split into 5-mers for the current analysis.

<https://www.nature.com/articles/s41598-018-21021-9>

## Deconvolution: Generation of Protease Reference Profile

A protease reference profile was constructed based on DESeq2-normalized
counts of cleaved peptides across the 7 conditions. Counts data from
neutrophil proteases, phage controls and ADAMTS13/thrombin were merged
and subsequently normalized via DESeq2.

``` r
#Reference Profile 1
ADAMTS13_thrombin <- read.csv("ADAMTS13_thrombin_5mer_data.csv", header=TRUE,row.names=1)
ADAMTS13_thrombin$peptide <- rownames(ADAMTS13_thrombin)
mergedData <- protease_counts
mergedData$peptide <- rownames(mergedData)
mergedData <- merge(mergedData, ADAMTS13_thrombin, by = "peptide", all = TRUE)
mergedData[is.na(mergedData)] <- 0
rownames(mergedData) <- mergedData$peptide
mergedData <- mergedData[, -1]

ADAMTS13_thrombin_conditions = c("cathepsin G","elastase","hPR3","no protease", "ADAMTS13", "thrombin")
ADAMTS13_thrombin_replicates = c("C1", "C2", "C3", "E1", "E2", "E3", "H1", "H2", "H3", "V1", "V2", "V3", "N1", "N2", "N3", "AD1", "AD2", "AD3", "AD4", "T1", "T2", "T3", "T4")
ADAMTS13_thrombin_col = c("cathepsin G","cathepsin G","cathepsin G", "elastase","elastase","elastase","hPR3","hPR3","hPR3","no protease", "no protease", "no protease", "FLAG","FLAG","FLAG", "ADAMTS13", "ADAMTS13", "ADAMTS13", "ADAMTS13", "thrombin", "thrombin", "thrombin", "thrombin")
ADAMTS13_thrombin_colname = c("cathepsin G 1","cathepsin G 2","cathepsin G 3","elastase 1","elastase 2","elastase 3","hPR3 1","hPR3 2","hPR3 3","no protease 1", "no protease 2", "no protease 3", "FLAG 1","FLAG 2","FLAG 3", "ADAMTS13 1", "ADAMTS13 2", "ADAMTS13 3", "ADAMTS13 4", "thrombin 1", "thrombin 2", "thrombin 3", "thrombin 4")
control = "FLAG"

data_NS <- subset(mergedData, grepl("X", row.names(mergedData))==FALSE)
countData <- as.matrix(data_NS)
colData <- data.frame(condition=factor(ADAMTS13_thrombin_col))
dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds$condition <- relevel(dds$condition, ref = "FLAG")
keep <- rowSums(counts(dds)) >= 10
dds_ADAMTS13_thrombin <- dds[keep,]
dds_ADAMTS13_thrombin=DESeq(dds_ADAMTS13_thrombin)
saveRDS(dds_ADAMTS13_thrombin, file = paste0("ADAMTS13_thrombin_dds_object.rds"))
```

## Deconvolution: Protease Counts Data Preparation

Following DESeq2 Counts Normalization, input data for CibersortX and
EPIC was prepared and formatted. Due to memory limitations on processing
by EPIC/CibersortX, a row sum filter of 140 was applied.

``` r
ADAMTS13_thrombin_matrix <- as.data.frame(counts(dds_ADAMTS13_thrombin, normalized = TRUE, replaced = FALSE))
row_names <- as.list(rownames(ADAMTS13_thrombin_matrix))
ADAMTS13_thrombin_matrix$condition <- row_names
ADAMTS13_thrombin_matrix <- ADAMTS13_thrombin_matrix[, c(ncol(ADAMTS13_thrombin_matrix), 1:(ncol(ADAMTS13_thrombin_matrix)-1))]
colnames(ADAMTS13_thrombin_matrix) <- c("peptide", as.list(ADAMTS13_thrombin_col))

ADAMTS13_counts <- as.data.frame(counts(dds_ADAMTS13_thrombin, replaced = FALSE))
ADAMTS13_counts <- ADAMTS13_counts[apply(ADAMTS13_counts, 1, function(row) sum(row) < 140), ]
ADAMTS13_thrombin_matrix <- ADAMTS13_thrombin_matrix[!rownames(ADAMTS13_thrombin_matrix) %in% rownames(ADAMTS13_counts), ]
ADAMTS13_thrombin_matrix <- as.matrix(ADAMTS13_thrombin_matrix)
write.table(ADAMTS13_thrombin_matrix, file = paste0("filtered_protease_ref_ADAMTS13_thrombin.txt"), sep = "\t", row.names = FALSE)
```

### Deconvolution: Mixture Counts Data Preparation

The bulk mixture dataset, which consists of the activated neutrophil
supernatant across different conditions, was generated based on
normalized cleaved peptide counts. Due to memory limitations on
processing by EPIC/CibersortX, a row sum filter of 30 was applied.

``` r
#Bulk Mixture Profile
mixture_matrix <- as.data.frame(counts(dds_mixture, normalized = TRUE, replaced = FALSE))
row_names <- as.list(rownames(mixture_matrix))
mixture_matrix$condition <- row_names
mixture_matrix <- mixture_matrix[, c(ncol(mixture_matrix), 1:(ncol(mixture_matrix)-1))]
colnames(mixture_matrix) <- c("peptide", as.list(mixture_col))
mixture_matrix <- as.matrix(mixture_matrix)
write.table(mixture_matrix, file = paste0("supernatant_mixture.txt"), sep = "\t", row.names = FALSE)

lowcounts <- as.data.frame(counts(dds_mixture, replaced = FALSE))
lowcounts <- lowcounts[apply(lowcounts, 1, function(row) sum(row) < 30), ]
mixture_matrix_filtered <- mixture_matrix[!rownames(mixture_matrix) %in% rownames(lowcounts), ]
mixture_matrix_filtered <- as.matrix(mixture_matrix_filtered)
write.table(mixture_matrix_filtered, file = paste0("filtered_supernatant_mixture.txt"), sep = "\t", row.names = FALSE)
```

## Deconvolution of the Neutrophil Supernatant Mixture

## Deconvolution Results from CibersortX

For CibersortX, a protease signature matrix was first created based on
the reference profile, independently from mixture data. This matrix
extracted signature peptides that characterize the presence of tested
proteases/conditions.

### Parameters for Creation of the Protease Signature Matrix: (‘Signature Matrix Creation’)

Default (kappa = 999, q-value = 0.01, sampling = 0.5, min expression =
1)

Quantile Normalization enabled

Min-Max Barcode Peptides per Condition: adjusted for experiment
complexity (2000-3000)

This signature matrix was applied towards the bulk mixture to analyze
protease composition within each neutrophil supernatant.

### Parameters for Mixture Deconvolution: (’Impute Cell Fractions, Absolute Mode)

Default (batch correction disabled, quantile normalization disabled,
permutations = 1)

<div class="figure" style="text-align: center">

<img src="ADAMTHR_cibersortref.png" alt="CibersortX MEROPs Results" width="30%" /><img src="ADAMTHR_cibersort1.png" alt="CibersortX MEROPs Results" width="30%" /><img src="ADAMTHR_cibersort2.png" alt="CibersortX MEROPs Results" width="30%" />
<p class="caption">
CibersortX MEROPs Results
</p>

</div>

Total Peptides included in the Protease Signature Matrix: 11221

## Deconvolution Results from EPIC

In contrast to CibersortX, EPIC does not require the prior generation of
a specific signature matrix. The protease reference profile, of
normalized cleaved peptides, was used directly for deconvolution
analysis onto the mixture reference profile.

A signature peptide list for EPIC’s analyses was obtained from the
intersection between protease and mixture peptides (138224 peptides).

``` r
ADAMTS13_thrombin_common_peptides <- intersect(rownames(ADAMTS13_thrombin_matrix), rownames(mixture_matrix_filtered))
write.table(ADAMTS13_thrombin_common_peptides, file="ADAMTS13_thrombin_signature_peptides.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

supernatant <- read.table(file = "filtered_supernatant_mixture.txt", header=TRUE,row.names=1)
reference_ADAMThrombin <- read.table(file = "filtered_protease_ref_ADAMTS13_thrombin.txt", header=TRUE,row.names=1)
immunedeconv::deconvolute_epic_custom(
  supernatant,
  reference_ADAMThrombin,
  ADAMTS13_thrombin_common_peptides,
  genes_var = NULL,
  mrna_quantities = NULL)
```

<div class="figure" style="text-align: center">

<img src="ADAMTHR_epic1.png" alt="EPIC MEROPs results" width="20%" height="20%" /><img src="ADAMTHR_epic2.png" alt="EPIC MEROPs results" width="20%" height="20%" /><img src="ADAMTHR_epic3.png" alt="EPIC MEROPs results" width="20%" height="20%" />
<p class="caption">
EPIC MEROPs results
</p>

</div>

# Data from CibersortX Results

CibersortX estimates of proteolytic activity levels across mixture
conditions were visualized with barplots and heatmaps. Relative
percentages depicted activity as part of 100%, while absolute scores
depicted unstandardized levels of activity/presence.

## CibersortX: Generation of Figures for Percentages

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-49-1.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-49-2.png" width="40%" />

## CibersortX: Generation of Figures for Absolute Scores

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-50-1.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-50-2.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-50-3.png" width="40%" />

### Statistical Comparisons

Statistical comparisons using a Kruskal-Wallis test for differences and
2-sample t-test between conditions was performed for the absolute
scores. Results were visualized with boxplots.

``` r
conditions <- c("cathepsin_G", "elastase", "hPR3", "no_protease", "FLAG_peptide", "ADAMTS13", "thrombin")
names = c("cathepsin G", "elastase",  "hPR3", "no protease", "FLAG peptide", "ADAMTS13", "thrombin")
condition_comparisons <- list( c("no inhibitor", "AEBSF"), c("no inhibitor", "AEBSF+EDTA"), c("no inhibitor", "FLAG control") )
graph_scale = c(0.5, 5, 2, 4, 6, 1, 1) #scale graph to fit y-axis

for(i in 1:length(conditions)){
  box_plot <- ggboxplot(cibersort_absolute_ADAMThrombin, x = "Condition", y = conditions[i], add = "jitter") + 
  labs(x = "supernatant condition", y = paste0("absolute scores for ", names[i])) +
  stat_compare_means(comparisons = condition_comparisons, method = "t.test") + 
  stat_compare_means(label.y = graph_scale[i])    
  show(box_plot)
}
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-1.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-2.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-3.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-4.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-5.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-6.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-51-7.png" width="40%" />

# Data from EPIC Results

EPIC estimates of proteolytic activity levels across mixture conditions
were visualized with barplots and heatmaps. Relative percentages
depicted activity as part of 100%. Absolute scores were not reported by
EPIC.

## EPIC: Generation of Figures for Percentages

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-52-1.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-52-2.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-52-3.png" width="40%" />

### Statistical Comparisons

Statistical comparisons using a Kruskal-Wallis test for differences and
2-sample t-test between conditions was performed for the percentages.

``` r
conditions <- c("cathepsin_G", "elastase", "hPR3", "no_protease", "FLAG_peptide", "ADAMTS13", "thrombin", "other")
names = c("cathepsin G", "elastase","hPR3","no protease", "FLAG peptide", "ADAMTS13", "thrombin", "other")
condition_comparisons <- list( c("no inhibitor", "AEBSF"), c("no inhibitor", "AEBSF+EDTA"), c("no inhibitor", "FLAG control") )
graph_scale = c(5, 90, 60, 70, 80, 5, 5, 5) #scale graph to fit y-axis

for(i in 1:length(conditions)){
  box_plot <- ggboxplot(epic_percentage_ADAMThrombin, x = "Condition", y = conditions[i], add = "jitter") + 
  labs(x = "supernatant condition", y = paste0("absolute scores for ", names[i])) +
  stat_compare_means(comparisons = condition_comparisons, method = "t.test") + 
  stat_compare_means(label.y = graph_scale[i])    
  show(box_plot)
}
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-1.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-2.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-3.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-4.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-5.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-6.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-7.png" width="40%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-53-8.png" width="40%" />

# Repeated Deconvolution Analyses with different protease references

To test the robustness of the algorithms and obtained deconvolution
estimates, CibersortX and EPIC were repeated with different protease
reference profiles. Due to concerns specified in the manuscript, a trial
was conducted with hPR3 excluded.

Reference Profile 2) 3 neutrophil serine proteases, 2 phage controls
Reference Profile 3) cathepsin G, elastase, 2 phage controls (without
hPR3)

## Deconvolution: Protease Reference Matrix

A protease reference profile of the 3 neutrophil serine proteases and 2
phage controls was constructed. This was used as input for CibersortX to
construct the protease signature matrix.

### CibersortX Parameters for Protease Signature Matrix:

Default (kappa = 999, q-value = 0.01, sampling = 0.5, min expression =
1)

Quantile Normalization enabled

Min-Max Barcode Peptides per Condition: adjusted for experiment
complexity (500-1000 for all proteases, 2000-3000 without hPR3)

``` r
#Reference Profile 2
protease_matrix <- as.data.frame(counts(dds_protease, normalized = TRUE, replaced = FALSE))
row_names <- as.list(rownames(protease_matrix))
protease_matrix$condition <- row_names
protease_matrix <- protease_matrix[, c(ncol(protease_matrix), 1:(ncol(protease_matrix)-1))]
colnames(protease_matrix) <- c("peptide", as.list(protease_col))
protease_matrix <- as.matrix(protease_matrix)
write.table(protease_matrix, file = paste0("protease_ref.txt"), sep = "\t", row.names = FALSE)
```

``` r
#Reference Profile 3
protease_matrix_nohPR3 <- protease_matrix[, !grepl("^hPR3", colnames(protease_matrix))]
protease_matrix_nohPR3 <- as.matrix(protease_matrix_nohPR3)
write.table(protease_matrix_nohPR3, file = paste0("protease_ref_hPR3_excluded.txt"), sep = "\t", row.names = FALSE)
```

Total Peptides included in Signature Matrix: 6259 (All proteases), 6166
(hPR3 excluded)

<div class="figure" style="text-align: center">

<img src="referencenohPR3.png" alt="protease Signature Matrices without (left) and with (right) the inclusion of hPR3" width="30%" height="20%" /><img src="referencewithhPR3.png" alt="protease Signature Matrices without (left) and with (right) the inclusion of hPR3" width="30%" height="20%" />
<p class="caption">
protease Signature Matrices without (left) and with (right) the
inclusion of hPR3
</p>

</div>

### CibersortX Parameters for Mixture Deconvolution:

Protease composition analysis was based on default parameters: Default
(Batch Correction disabled, Quantile Normalization disabled,
Permutations = 1)

<div class="figure" style="text-align: center">

<img src="hpr3cibersortplot.png" alt="protease Percentages (left) and absolute scores (right) for all proteases" width="30%" height="50%" /><img src="mixturetable_withhPR3.png" alt="protease Percentages (left) and absolute scores (right) for all proteases" width="30%" height="50%" />
<p class="caption">
protease Percentages (left) and absolute scores (right) for all
proteases
</p>

</div>

<div class="figure" style="text-align: center">

<img src="mixturetable_withouthPR3.png" alt="protease Percentages (left) and absolute scores (right) for proteases with hPR3 excluded" width="30%" height="50%" /><img src="nohpr3cibersortplot.png" alt="protease Percentages (left) and absolute scores (right) for proteases with hPR3 excluded" width="30%" height="50%" />
<p class="caption">
protease Percentages (left) and absolute scores (right) for proteases
with hPR3 excluded
</p>

</div>

### EPIC Parameters for Protease Reference and Deconvolution:

A signature peptide list for EPIC’s analyses was obtained from the
intersection between protease and mixture peptides (203974 peptides).

Due to memory limitations on processing by EPIC , a row sum filter of 30
was applied.

``` r
lowcounts <- as.data.frame(counts(dds_protease, replaced = FALSE))
lowcounts <- lowcounts[apply(lowcounts, 1, function(row) sum(row) < 30), ]
matrix_filtered <- protease_matrix[!rownames(protease_matrix) %in% rownames(lowcounts), ]
matrix_filtered <- as.matrix(matrix_filtered)

matrix_filtered_nohPR3 <- matrix_filtered[, !grepl("^hPR3", colnames(matrix_filtered))]
matrix_filtered_nohPR3 <- as.matrix(matrix_filtered_nohPR3)
write.table(matrix_filtered, file = paste0("filtered_protease_ref.txt"), sep = "\t", row.names = FALSE)
write.table(matrix_filtered_nohPR3, file = paste0("filtered_protease_ref_nohPR3.txt"), sep = "\t", row.names = FALSE)

common_peptides <- intersect(rownames(matrix_filtered), rownames(mixture_matrix_filtered))
write.table(common_peptides, file="signature_peptides.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

<div class="figure" style="text-align: center">

<img src="hpr3EPIC1.png" alt="EPIC Results with all proteases" width="30%" height="20%" /><img src="hpr3EPIC2.png" alt="EPIC Results with all proteases" width="30%" height="20%" /><img src="hpr3EPIC3.png" alt="EPIC Results with all proteases" width="30%" height="20%" />
<p class="caption">
EPIC Results with all proteases
</p>

</div>

<div class="figure" style="text-align: center">

<img src="nohpr3EPIC1.png" alt="EPIC Results without hPR3" width="30%" height="20%" /><img src="nohpr3EPIC2.png" alt="EPIC Results without hPR3" width="30%" height="20%" /><img src="nohpr3EPIC3.png" alt="EPIC Results without hPR3" width="30%" height="20%" />
<p class="caption">
EPIC Results without hPR3
</p>

</div>

# Data from CibersortX Results: all neutrophil serine proteases (profile 2)

To examine the estimated protease % composition and absolute presence
score across mixture conditions from CibersortX, barplots and heatmaps
were generated.

## CibersortX: Generation of Figures for Percentages

``` r
percentage_data <- cibersort_percentage %>%
  mutate(Condition = paste0(Condition, " ", (row_number()-1)%%3 + 1))
rownames(percentage_data) <- percentage_data[, 1]
heatmap_data <- percentage_data[, -1]
pheatmap(heatmap_data, display_numbers = TRUE)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-63-1.png" width="30%" />

``` r
barplot_data <- reshape2::melt(percentage_data, id.vars = "Condition")
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = value, fill = variable)) + geom_bar(stat = "identity") + labs(x = "supernatant condition", y = "relative percentages", fill = "protease") + scale_fill_viridis_d()+ theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-63-2.png" width="30%" />

## CibersortX: Generation of Figures for Absolute Scores

``` r
absolute_data <- cibersort_absolute %>%
  mutate(Condition = paste0(Condition, " ", (row_number()-1)%%3 + 1))
barplot_data <- reshape2::melt(absolute_data, id.vars = "Condition")
rownames(absolute_data) <- absolute_data[, 1]
pheatmap(absolute_data[,-1], display_numbers = TRUE)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-64-1.png" width="30%" />

``` r
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = value, fill = variable)) + geom_bar(stat = "identity") + labs(x = "supernatant condition", y = "absolute scores", fill = "protease") + scale_fill_viridis_d()+ theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-64-2.png" width="30%" />

``` r
barplot_data <- cibersort_absolute %>%
  pivot_longer(-Condition, names_to = "protease", values_to = "Absolute_Score")
barplot_data <- barplot_data %>% 
mutate(protease = factor(protease, levels = c("cathepsin_G", "elastase", "hPR3", "no_protease", "FLAG_peptide")))
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = Absolute_Score, fill = protease)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(x = "supernatant condition", y = "absolute scores", fill = "protease") +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_beeswarm(data = barplot_data, aes(color = protease), size = 2, alpha = 0.5, position = position_dodge(width = 0.9))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-64-3.png" width="30%" />

### Statistical Comparisons

``` r
conditions <- c("cathepsin_G", "elastase", "hPR3", "no_protease", "FLAG_peptide")
names = c("cathepsin G", "elastase", "hPR3", "no protease", "FLAG peptide")
condition_comparisons <- list( c("no inhibitor", "AEBSF"), c("no inhibitor", "AEBSF+EDTA"), c("no inhibitor", "FLAG control") )
graph_scale = c(5, 45, 10, 25, 35) #scale graph to fit y-axis

for(i in 1:length(conditions)){
  box_plot <- ggboxplot(cibersort_absolute, x = "Condition", y = conditions[i], add = "jitter") + 
  labs(x = "supernatant condition", y = paste0("absolute scores for ", names[i])) +
  stat_compare_means(comparisons = condition_comparisons, method = "t.test") + 
  stat_compare_means(label.y = graph_scale[i])    
  show(box_plot)
}
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-65-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-65-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-65-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-65-4.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-65-5.png" width="30%" />

# Data from EPIC Results: all neutrophil serine proteases (profile 2)

To examine the estimated protease % composition and absolute presence
score across mixture conditions from CibersortX, barplots and heatmaps
were generated.

## EPIC: Generation of Figures for Percentages

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-67-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-67-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-67-3.png" width="30%" />

### Statistical Comparisons

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-4.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-5.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-68-6.png" width="30%" />

# Data from CibersortX Results: cathepsin G and elastase (profile 3)

Analyses from previously were repeated.

## CibersortX: Generation of Figures for Percentages

``` r
percentage_data <- cibersort_percentage_nohPR3 %>%
  mutate(Condition = paste0(Condition, " ", (row_number()-1)%%3 + 1))
rownames(percentage_data) <- percentage_data[, 1]
heatmap_data <- percentage_data[, -1]
pheatmap(heatmap_data, display_numbers = TRUE)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-69-1.png" width="30%" />

``` r
barplot_data <- reshape2::melt(percentage_data, id.vars = "Condition")
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = value, fill = variable)) + geom_bar(stat = "identity") + labs(x = "supernatant condition", y = "relative percentages", fill = "protease") + scale_fill_viridis_d()+ theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-69-2.png" width="30%" />

## CibersortX: Generation of Figures for Absolute Scores

``` r
absolute_data <- cibersort_absolute_nohPR3 %>%
  mutate(Condition = paste0(Condition, " ", (row_number()-1)%%3 + 1))
barplot_data <- reshape2::melt(absolute_data, id.vars = "Condition")
rownames(absolute_data) <- absolute_data[, 1]
pheatmap(absolute_data[,-1], display_numbers = TRUE)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-70-1.png" width="30%" />

``` r
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = value, fill = variable)) + geom_bar(stat = "identity") + labs(x = "supernatant condition", y = "absolute scores", fill = "protease") + scale_fill_viridis_d()+ theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-70-2.png" width="30%" />

``` r
barplot_data <- cibersort_absolute_nohPR3 %>%
  pivot_longer(-Condition, names_to = "protease", values_to = "Absolute_Score")
barplot_data <- barplot_data %>% 
  mutate(protease = factor(protease, levels = c("cathepsin_G", "elastase", "no_protease", "FLAG_peptide")))
ggplot(barplot_data, aes(x = fct_inorder(Condition), y = Absolute_Score, fill = protease)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9)) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(x = "supernatant condition", y = "absolute scores", fill = "protease") +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_beeswarm(data = barplot_data, aes(color = protease), size = 2, alpha = 0.5, position = position_dodge(width = 0.9))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-70-3.png" width="30%" />

### Statistical Comparisons

``` r
conditions <- c("cathepsin_G", "elastase", "no_protease", "FLAG_peptide")
names = c("cathepsin G", "elastase","no protease", "FLAG peptide")
condition_comparisons <- list( c("no inhibitor", "AEBSF"), c("no inhibitor", "AEBSF+EDTA"), c("no inhibitor", "FLAG control") )
graph_scale = c(10, 80, 25, 45) #scale graph to fit y-axis

for(i in 1:length(conditions)){
  box_plot <-ggboxplot(cibersort_absolute_nohPR3, x = "Condition", y = conditions[i], add = "jitter") + 
  labs(x = "supernatant condition", y = paste0("absolute scores for ", names[i])) +
  stat_compare_means(comparisons = condition_comparisons, method = "t.test") + 
  stat_compare_means(label.y = graph_scale[i])    
  show(box_plot)
}
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-71-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-71-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-71-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-71-4.png" width="30%" />

# Data from EPIC Results: cathepsin G and elastase (profile 3)

To examine the estimated protease % composition and absolute presence
score across mixture conditions from CibersortX, barplots and heatmaps
were generated.

## EPIC: Generation of Figures for Percentages

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-72-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-72-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-72-3.png" width="30%" />

### Statistical Comparisons

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-73-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-73-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-73-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-73-4.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-73-5.png" width="30%" />

# Comparison to MEROPs Database

The MEROPs database consists of known substrates for individual
proteases, like neutrophil serine proteases, based on prior studies and
screens. To validate the current dataset towards known substrates of
proteases, its cleaved peptides were compared to those of MEROPs.

MEROPs substrates (8-mers) were obtained and converted into 5-mer or
4-mer space, which were subsequently compared with the significantly
cleaved peptides of cathepsin G, elastase, and hPR3.

``` r
convert_3_to_1 <- function(three_letter_code) {
  codes <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V"
  )
  if(is.na(codes[three_letter_code])){
    return("-")
  } else{
    return(codes[three_letter_code])
  }
}
```

## Splitting of MEROPs 8-mers into 5-mers/4-mers

A MEROPs database list of substrates was obtained for each of the 3
neutrophil serine proteases: cathepsin G, neutrophil elastase, and human
proteinase 3. Tables were derived and 5-mers/4-mers were compiled for
each 8-mer substrate. Only 5-mers/4-mers with natural amino acids were
included for analysis.

``` r
meropsList <- list()
meropsList_4mer <- list()
combined <- append(proteases, mixtures)

for(i in 1:3){
  meropsData <- read.xlsx(paste0(combined[i], "_MEROPs.xlsx"))
  aminoacidMerops <- meropsData[7:14]
  colnames(aminoacidMerops) <- c("P4",  "P3",  "P2",  "P1", "P1x", "P2x", "P3x","P4x")
  aminoacidMerops$P4 <- sapply(aminoacidMerops$P4, convert_3_to_1) 
  aminoacidMerops$P3 <- sapply(aminoacidMerops$P3, convert_3_to_1) 
  aminoacidMerops$P2 <- sapply(aminoacidMerops$P2, convert_3_to_1) 
  aminoacidMerops$P1 <- sapply(aminoacidMerops$P1, convert_3_to_1) 
  aminoacidMerops$P1x <- sapply(aminoacidMerops$P1x, convert_3_to_1)
  aminoacidMerops$P2x <- sapply(aminoacidMerops$P2x, convert_3_to_1) 
  aminoacidMerops$P3x <- sapply(aminoacidMerops$P3x, convert_3_to_1) 
  aminoacidMerops$P4x <- sapply(aminoacidMerops$P4x, convert_3_to_1)
  write.table(aminoacidMerops, file = paste0(combined[i], "_merops_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  short_mer <- c()
  for(j in 1:length(aminoacidMerops[,1])){
    row <- aminoacidMerops[j,]
    for(k in 1:5){
      mer <- apply(matrix(row, nrow = 4), 2, paste, collapse = "")[1]
      row <- row[,-1]
      if(!grepl("-", mer)){
       short_mer <- c(short_mer, mer)
      }
    }
  }
  write.table(as.data.frame(short_mer), file = paste0(combined[i], "_merops_4mer_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
   meropsList_4mer <- append(meropsList_4mer, list(short_mer))
  names(meropsList_4mer)[length(meropsList_4mer)] <- combined[i]
  short_mer <- c()
  for(j in 1:length(aminoacidMerops[,1])){
    row <- aminoacidMerops[j,]
    for(k in 1:4){
      mer <- apply(matrix(row, nrow = 5), 2, paste, collapse = "")[1]
      row <- row[,-1]
      if(!grepl("-", mer)){
       short_mer <- c(short_mer, mer)
      }
    }
  }
  write.table(as.data.frame(short_mer), file = paste0(combined[i], "_merops_5mer_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  meropsList <- append(meropsList, list(short_mer))
  names(meropsList)[length(meropsList)] <- combined[i]
}
```

MEROPs substrates for each protease were visualized using ggseqlogo.

``` r
ggseqlogo(meropsList, ncol = 5)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-76-1.png" width="40%" />

``` r
ggseqlogo(meropsList_4mer, ncol = 5)
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-76-2.png" width="40%" />

## Generation of 4-mer Significantly Cleaved Peptides

Significantly cleaved peptides from the phage dataset were shifted into
4-mer space, for subsequent comparison to MEROPs 4-mers.

``` r
peptide_4mers <- list()
for(i in 1:3){
  seqs <- peptideList[[i]]
  temp <- c()
  for(j in 1:length(seqs)){
    temp <- c(temp, substr(seqs[j], 1, 4),substr(seqs[j], 2, 5))
  }
  peptide_4mers <- append(peptide_4mers, list(temp))
  names(peptide_4mers)[length(peptide_4mers)] <- combined[i]
}
```

## Venn Diagram Comparisons between 5-mers

To identify similarities and differences in motifs, 5-mer peptides of
MEROPs substrates were compared with phage cleaved peptides using
ggVennDiagram.

``` r
venn_list <- list(`cathepsin G MEROPs 5mer` = meropsList[[1]], `elastase MEROPs 5mer`=meropsList[[2]], `hPR3 MEROPs 5mer`=meropsList[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-78-1.png" width="40%" />

``` r
venn_list <-  list(`experiment cathepsin G` = peptideList[[1]], `cathepsin G MEROPs 5mer` = meropsList[[1]],`elastase MEROPs 5mer`=meropsList[[2]], `hPR3 MEROPs 5mer`=meropsList[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-78-2.png" width="40%" />

``` r
venn_list <-  list(`experiment elastase`= peptideList[[2]], `cathepsin G MEROPs 5mer` = meropsList[[1]],`elastase MEROPs 5mer`=meropsList[[2]], `hPR3 MEROPs 5mer`=meropsList[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-78-3.png" width="40%" />

``` r
venn_list <-  list(`experiment hPR3` = peptideList[[3]], `cathepsin G MEROPs 5mer` = meropsList[[1]],`elastase MEROPs 5mer`=meropsList[[2]], `hPR3 MEROPs 5mer`=meropsList[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-78-4.png" width="40%" />

## Venn Diagram Comparisons between 4-mers

To identify similarities and differences in motifs, 4-mer peptides of
MEROPs substrates were compared with phage cleaved peptides using
ggVennDiagram.

``` r
venn_list <- list(`cathepsin G MEROPs 4mer` = meropsList_4mer[[1]],`elastase MEROPs 4mer`=meropsList_4mer[[2]], `hPR3 MEROPs 4mer`=meropsList_4mer[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-79-1.png" width="40%" />

``` r
venn_list <-  list(`experimental cathepsin G` = peptide_4mers[[1]], `cathepsin G MEROPs 4mer` = meropsList_4mer[[1]],`elastase MEROPs 4mer`=meropsList_4mer[[2]], `hPR3 MEROPs 4mer`=meropsList_4mer[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-79-2.png" width="40%" />

``` r
venn_list <-  list(`experimental elastase` = peptide_4mers[[2]], `cathepsin G MEROPs 4mer` = meropsList_4mer[[1]],`elastase MEROPs 4mer`=meropsList_4mer[[2]], `hPR3 MEROPs 4mer`=meropsList_4mer[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-79-3.png" width="40%" />

``` r
venn_list <-  list(`experimental hPR3` = peptide_4mers[[3]], `cathepsin G MEROPs 4mer` = meropsList_4mer[[1]],`elastase MEROPs 4mer`=meropsList_4mer[[2]], `hPR3 MEROPs 4mer`=meropsList_4mer[[3]])
ggVennDiagram(venn_list) + scale_x_continuous(expand = expansion(mult = .2))
```

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-79-4.png" width="40%" />

# Deconvolution Analysis using MEROPs-Based Protease Reference Matrices

To assess the applicability of existing MEROPs substrate data toward
previous deconvolution analyses of complex biological mixtures, MEROPs
5-mer frequencies were applied as a signature protease matrix for
CibersortX and EPIC. As such, instead of the phage display’s protease
cleaved peptides, MEROPs-derived peptides were used. Other conditions
and data remained consistent.

The MEROPs analysis was repeated with and without the FLAG and no
protease controls generated from the phage experiments.

Reference Profile 4) 3 neutrophil serine proteases from MEROPs
substrates, 2 phage controls Reference Profile 5) 3 neutrophil serine
proteases from MEROPs substrates

## Generation of the MEROPs reference matrix

``` r
cathepsin_G <- as.data.frame(meropsList[1])
cathepsin_G <- table(unlist(cathepsin_G))
cathepsin_G <- as.data.frame(cathepsin_G)
colnames(cathepsin_G) <- c("peptide", "frequency")

elastase <- as.data.frame(meropsList[2])
elastase <- table(unlist(elastase))
elastase <- as.data.frame(elastase)
colnames(elastase) <- c("peptide", "frequency")

hPR3 <- as.data.frame(meropsList[3])
hPR3 <- table(unlist(hPR3))
hPR3 <- as.data.frame(hPR3)
colnames(hPR3) <- c("peptide", "frequency")

MEROPs_df <- merge(merge(cathepsin_G, elastase, by = "peptide", all = TRUE), hPR3, by = "peptide", all = TRUE)
phage_control_data <- protease_counts[c(10, 11, 12, 13, 14, 15)]
phage_control_data$peptide <- rownames(phage_control_data)
phage_control_data <- phage_control_data[, c(ncol(phage_control_data), 1:(ncol(phage_control_data)-1))]
MEROPs_and_control <- merge(MEROPs_df, phage_control_data, by = "peptide", all = TRUE)
MEROPs_df[is.na(MEROPs_df)] <- 0
MEROPs_and_control[is.na(MEROPs_and_control)] <- 0
```

## Normalization and Processing for Deconvolution Input

The following describes results from the previous analyses, with and
without the phage controls.

Normalized counts data was generated for CibersortX and EPIC input,
scaled up by 2000. Data for EPIC was rounded to 4 decimal places for
processing. Common peptides between the supernatant and MEROPs reference
matrix were selected as signature peptides.

``` r
#Reference Profile 5
MEROPs_df$peptide <- as.character(MEROPs_df$peptide)
MEROPs_df <- MEROPs_df[order(MEROPs_df$peptide), ]
normalized_data <- scale(MEROPs_df[,-1], center = FALSE, scale = colSums(MEROPs_df[, -1]))
normalized_data <- as.data.frame(normalized_data)*2000
normalized_df <- cbind(MEROPs_df[, 1, drop = FALSE], as.data.frame(normalized_data))
repeated_df <- cbind(normalized_df, normalized_df[, 2], normalized_df[, 3], normalized_df[, 4], normalized_df[, 2], normalized_df[, 3], normalized_df[, 4])
colnames(repeated_df) <- c("peptide", "cathepsin_G", "elastase", "hPR3", "cathepsin_G", "elastase", "hPR3", "cathepsin_G", "elastase", "hPR3")
write.table(repeated_df, file = paste0("protease_ref_MEROPs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

MEROPs_common_peptides <- intersect(repeated_df$peptide, rownames(mixture_matrix_filtered))
write.table(MEROPs_common_peptides, file="MEROPs_signature_peptides.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Reference Profile 4
MEROPs_and_control$peptide <- as.character(MEROPs_and_control$peptide)
MEROPs_and_control <- MEROPs_and_control[order(MEROPs_and_control$peptide), ]
normalized_data <- scale(MEROPs_and_control[,-1], center = FALSE, scale = colSums(MEROPs_and_control[, -1]))
normalized_data <- as.data.frame(normalized_data)*2000
normalized_df <- cbind(MEROPs_and_control[, 1, drop = FALSE], as.data.frame(normalized_data))
repeated_df <- cbind(normalized_df, normalized_df[, 2], normalized_df[, 3], normalized_df[, 4], normalized_df[, 2], normalized_df[, 3], normalized_df[, 4])
colnames(repeated_df) <- c("peptide", "cathepsin_G", "elastase", "hPR3", "no_protease", "no_protease", "no_protease", "FLAG_peptide", "FLAG_peptide", "FLAG_peptide", "cathepsin_G", "elastase", "hPR3", "cathepsin_G", "elastase", "hPR3")
write.table(repeated_df, file = paste0("protease_ref_MEROPs_control.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

round_df <- function(x, digits) {
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x}
shrunk_df <- round_df(repeated_df, 4)
write.table(shrunk_df, file = paste0("shrunk_protease_ref_MEROPs_control.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

MEROPs_common_peptides_control <- intersect(shrunk_df$peptide, rownames(mixture_matrix_filtered))
write.table(MEROPs_common_peptides_control, file="MEROPs_signature_peptides_control.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

## Deconvolution of Mixtures with MEROPs-based protease references and phage controls (profile 4)

### Deconvolution Results from CibersortX

Protease Signature Matrix Parameters:

Default (kappa = 999, q-value = 0.01, sampling = 0.5, min expression =
1)

Quantile Normalization enabled

Min-Max Barcode Peptides per Condition: adjusted for experiment
complexity (500-1000 for all proteases with controls and all proteases
without controls)

Total Peptides included in Signature Matrix: 1041 (all proteases with
controls), 82 (all proteases without controls)

Deconvolution Fraction Parameters:

Default (batch correction disabled, quantile normalization disabled,
permutations = 1)

<div class="figure" style="text-align: center">

<img src="cibermerops_controlref.png" alt="CibersortX MEROPs results with phage control" width="30%" height="20%" /><img src="cibermerops_control.png" alt="CibersortX MEROPs results with phage control" width="30%" height="20%" /><img src="cibermerops_control2.png" alt="CibersortX MEROPs results with phage control" width="30%" height="20%" />
<p class="caption">
CibersortX MEROPs results with phage control
</p>

</div>

Generated barplots and heatmaps from results

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-84-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-84-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-84-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-84-4.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-84-5.png" width="30%" />

### Deconvolution Results from EPIC

EPIC signature peptides: 207119 (all proteases with controls), 283 (all
proteases without controls)

<div class="figure" style="text-align: center">

<img src="epicmerops_control1.png" alt="EPIC MEROPs results with phage control" width="30%" height="20%" /><img src="epicmerops_control2.png" alt="EPIC MEROPs results with phage control" width="30%" height="20%" /><img src="epicmerops_control3.png" alt="EPIC MEROPs results with phage control" width="30%" height="20%" />
<p class="caption">
EPIC MEROPs results with phage control
</p>

</div>

Generated barplots and heatmaps from results

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-86-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-86-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-86-3.png" width="30%" />

## Deconvolution of Mixtures with MEROPs-based protease references and no phage controls (profile 5)

### Deconvolution Results from CibersortX

<div class="figure" style="text-align: center">

<img src="cibermerops_nocontrolref.png" alt="CibersortX MEROPs results without phage control" width="30%" height="20%" /><img src="cibermerops_nocontrol.png" alt="CibersortX MEROPs results without phage control" width="30%" height="20%" /><img src="cibermerops_nocontrol2.png" alt="CibersortX MEROPs results without phage control" width="30%" height="20%" />
<p class="caption">
CibersortX MEROPs results without phage control
</p>

</div>

Generated barplots and heatmaps from results

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-89-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-89-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-89-3.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-89-4.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-89-5.png" width="30%" />

### Deconvolution Results from EPIC

<div class="figure" style="text-align: center">

<img src="epicmerops_nocontrol1.png" alt="EPIC MEROPs results without phage control" width="30%" height="20%" /><img src="epicmerops_nocontrol2.png" alt="EPIC MEROPs results without phage control" width="30%" height="20%" /><img src="epicmerops_nocontrol3.png" alt="EPIC MEROPs results without phage control" width="30%" height="20%" />
<p class="caption">
EPIC MEROPs results without phage control
</p>

</div>

Generated barplots and heatmaps from results

<img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-91-1.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-91-2.png" width="30%" /><img src="Protease-Analysis-Mid-June_files/figure-gfm/unnamed-chunk-91-3.png" width="30%" />
