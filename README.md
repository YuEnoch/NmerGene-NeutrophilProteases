# NNK-NeutrophilProteases
**High Throughput Sequencing Analysis Pipeline for analyzing 5 Amino Acid Mutagenesis Data**

This folder contains files involved in the high throughput sequencing analysis from the 5 Amino Acid Mutagenesis sequencing data. It identifies significantly enriched peptides, analyzes the amino/nucleic acid frequencies, and clusters the enriched peptides based on their properties (Principal Component Analysis). The pipeline is all outlined in order within SeqSpacer.py, which automatically calls the rest of the scripts as long as they are in the same folder. 

This program can accommodate any n-mer mutagenesis library at any location, for any number of experments/replicates. Please specify these details in parameter.txt, which has all the default settings. It analyzes one-end reads, which could be single or in parallel. 

**To Start (based on default parameters, testing data):**
1. Download Repository and Unzip
2. Run SeqSpacer.py (*Ensure correct R Version)
3. A Python Popup should appear and the Pipeline would automatically run using the Test Files
4. The Plots and Reults would all appear within the file and its folders

**Prerequisites:**
- Python 3
- R (with R Script Path setup) 
- Installation of R Packages (specified in DESEQprocessor.R, PCAprocessor.R)
- Installation of Weblogo: https://github.com/WebLogo/weblogo, go to CMD and enter pip install weblogo

**To Adjust for different Sequencing Experiments:**
1. Change Parameters in parameters.txt (# of treatments, # of n_plicates, experiment/treatment names, gzip, etc.)
   - Ensure that the formatting is the same
2. Change Seeds (Conserved Sequences) in clean_fastq.py, based on your sequencing

**It contains:**

- Information about raw fastq file and its details
- All in-house developed scripts for read mapping, filtering, and translation
- DESEQ2 analysis for significantly enriched peptides
- Principal Component Analysis to cluster significatly enriched peptides
- Weblogo Analysis for clusters

- Frequency Analysis for Amino Acids, Nucleic Acids, and Oberved/Expected Amino Acids
- Code to generate plots for the above analyses
- Test Files (a subset of the complete FASTQ sequencing results), for Neutrophil Proteases Cathepsin G (C), Elastase (E), and Human Proteinase 3 (H) with an unselected control (N). The experiments were done in duplicates.

The code adapts work from:
- Dr. Kart Tomberg https://github.com/tombergk/NNK_VWF73/
- Dr. Matt Holding


mod575.py is written by Dr. Kart Tomberg

files in PCA Analysis from Dr. Matt Holding
