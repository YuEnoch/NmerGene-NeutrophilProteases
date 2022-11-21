# NNK-NeutrophilProteases
**High Throughput Sequencing Analysis Pipeline for analyzing 5 Amino Acid Mutagenesis Data**

This folder contains files involved in the high throughput sequencing analysis from the 5 Amino Acid Mutagenesis sequencing data. It identifies significantly enriched peptides, analyzes the amino/nucleic acid frequencies, and clusters the enriched peptides based on their properties (Principal Component Analysis).

The pipeline is all outlined in order within MainFile.py. To start, copy all Python and R scripts into the same folder. 

It contains:

- Information about raw fastq file and its details
- All in-house developed scripts for read mapping, filtering, and translation
- DESEQ2 analysis for significantly enriched peptides
- Principal Component Analysis to cluster significatly enriched peptides
- Weblogo Analysis for clusters

- Frequency Analysis for Amino Acids, Nucleic Acids, and Oberved/Expected Amino Acids
- Code to generate plots for the above analyses
- Test Files (a subset of the complete FASTQ sequencing results), for Neutrophil Proteases Cathepsin G (C), Human Proteinase 3 (H), and Elastase (E), with an unselected control (N). The tests were done in duplicates.

The code adapts work from:
- Dr. Kart Tomberg
- Dr. Matt Holding


mod575.py is written by Dr. Kart Tomberg

files in PCA Analysis from Dr. Matt Holding
