# High Throughput Sequencing Analysis Pipeline for analyzing N-Mer Amino Acid Mutagenesis Data

This folder contains files involved in the high throughput sequencing analysis from the 5 Amino Acid Mutagenesis sequencing data. It identifies significantly enriched peptides, analyzes the amino/nucleic acid frequencies, and clusters the enriched peptides based on their properties (Principal Component Analysis). The pipeline is all outlined in order within NmerGene.py, which automatically calls the rest of the scripts as long as they are in the same folder. 

This program can accommodate any n-mer mutagenesis library at any location, for any number of experments/replicates. Please specify these details in parameter.txt, which has all the default settings. It analyzes one-end reads, which could be single or in parallel. 

**To Start (based on default parameters, testing data):**
1. Download Repository and Unzip Folder
2. Run NmerGene.py
3. A Python Popup should appear and the Pipeline would automatically run using the Test Files
4. The Plots and Reults would all appear within the file and its folders

**Prerequisites:**
* Python 3
* R (with R Script Path setup) 
* [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [MClust](https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html)
* [Weblogo](https://github.com/WebLogo/weblogo)
* Installation of Other R Packages for Visualization (specified in the R Scripts)


**To Adjust for different Sequencing Experiments:**
1. Change Parameters in parameters.txt (# of treatments, # of n_plicates, experiment/treatment names, gzip, etc.)
   - Ensure that the formatting is the same
2. Change Seeds (Conserved Sequences) in clean_fastq.py, based on your sequencing

## Contents

* [Set-Up] (#set-up)
* [Raw FastQ File Processing - Trimming, Filtering, Translation] (#raw-fastq-file-processing---trimming-filtering-translation)
* [DESEQ2 Analysis for Significantly Enriched Peptides using R] (#deseq2)
* [Principal Component Analysis to cluster Significatly Enriched Peptides] (#pca)
* [Weblogo Analysis for clusters] (#weblogo)

* [Amino Acid Frequency Analysis] (#aa-freq)
* [Nucleic Acid Frequency Analysis] (#na-freq)
* [Observed/Expected Amino Acid Frequency Analysis] (#obs-exp)
* [Test Files] (#test-files)

### Set-Up
To set up this program, place all the files in the same folder (as organized by default). All folder directories would be set-up automatically, with the program ready to run on the test set (run NmerGene.py). The parameters.txt contains all the editable parameters available. Please edit each while maintaining the same formatting and spacing. 

Experiment Criteria: Any number of experiments, with 1 control. The number of n-plicates should be the same for all (e.g. A1, A2, B1, B2, N1, N2, where N is the control). This program accepts one-end reads (single or parallel). If they are paired-end, please merge them prior.

To run the program, place the raw FastQ or Gzip files within the folder and edit the parameters.txt accordingly (below). Once everything is saved, run NmerGene.py

```
*** Treatment Details ***
Experiment Name = NeutrophilProteases
Experiments = cathepsinG, elastase, hpr3
Control Name = N

Treatments = 8
N-plicates (duplicate/triplicate/etc.) = 2
Treatment Names = C, E, H, N
All Treatment Names = C1, C2, E1, E2, H1, H2, N1, N2

*** Mutagenesis Sequencing ***
First Position = 55
Last Position = 69
N-mer = 5
Min Count for each peptide (for DESEQ2) = 6

***
Files in Gzip (no/yes) = no
Single/Parallel One-End Read = parallel

*** Principal Component Analysis ***
Principal Components Used = 3
Limit of Clusters = 10 
```


### Raw FastQ File Processing - Trimming, Filtering, Translation



(a subset of the complete FASTQ sequencing results), for Neutrophil Proteases Cathepsin G (C), Elastase (E), and Human Proteinase 3 (H) with an unselected control (N). The experiments were done in duplicates.

The code adapts work from:
- Dr. Kart Tomberg https://github.com/tombergk/NNK_VWF73/
- Dr. Matt Holding


mod575.py is written by Dr. Kart Tomberg

files in PCA Analysis from Dr. Matt Holding
