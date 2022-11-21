#YuEnoch 13-11-2022
#MainFile.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: Runs All the Analyses involved, takes in 5AA Mutagenesis Data to derive a merged table (ready for Deseq2), AA_frequencies, Nucleic Acid frequencies
#Assumptions:
# 1. Mutagenesis for 5 Amino Acids
# 2. Mutagenesis occurs between positions 55 and 70 (in clean_fastq.py and get_peptides.py)
# 3. Matched Seeds (the conserved sequences used), in clean_fastq.py
# 4. Quality Threshold (here, removes reads when any two consecutive quality scores are both below 28)
# 5. For each peptide, there must be a total count over 6 across all treatments

#Requirements:
# 1. All FastQ files are stored in the same folder
# 2. Single-End Reads conducted in Parallel (each treatment has 2 FastQ Files, separate and not complementary)

#NOTE: need to edit the following (based off your experiment)
#      currently, all parameters are set based off the Test Dataset
treatments = 8
n_plicates = 2 #duplicates, triplicates, etc.
allTreatments = ["C1", "C2", "E1", "E2", "H1", "H2", "N1", "N2"]
treatmentNames = ["C", "E", "H", "N"]
experiments = ["cathepsinG", "elastase", "hpr3"]
experimentName = "Neutrophils"
gzipChoice = str(0) #Are the FastQ files in Gzip? 0 = No, 1 = Yes

def fileNameFormat(num, parallel):     #FastQ File Name Format
    return allTreatments[num] + "_S" + str(num+1) + "_L00"+str(parallel)+"_R1_001.fastq"

#NOTE: for the treatments and experiments, they need to be edited in the R files as well
#      the directories need to be changed in R files
#      also, the R files require packages to be installed (details in R files)


import sys
import os
import subprocess
from subprocess import call

def find_files(fileName):              #Find File Location (if FastQ not in same folder), goes outside 2 Folders to look for file
    dir_path = os.path.normpath(os.path.normpath(os.getcwd() + os.sep + os.pardir) + os.sep + os.pardir)   #goes outside two folders
    for root, dirs, files in os.walk(dir_path):
        if fileName in files: 
            return root+'/'+str(fileName)    

#clean_fastq.py: takes in raw FastQ data (single-end read) and filters sequencing reads by Quality and Matched Seeds, to obtain good Sense sequences (and Junk data)
#getPeptides.py: takes in Sense sequences and translates the Mutagenesis Region to obtain a table of Peptide Frequencies

for i in range(len(allTreatments)):
    name = allTreatments[i]
    fileName1 = find_files(fileNameFormat(i, 1))
    fileName2 = find_files(fileNameFormat(i, 2))
    print(name, fileName1, fileName2)
    
    call(["python3", "clean_fastq.py", name, fileName1, fileName2, gzipChoice])
    print(name, "clean_fastq.py done")
    call(["python3", "getPeptides.py", name])
    print(name, "getPeptides.py done")


    call(["python3", "getAminoAcidCount.py", name])
    call(["python3", "getNucleicAcidCount.py", name])
    call(["python3", "ObservedExpectedPlot.py", name])
    


#Generates Observed/Expected, Nucleotide, Amino Acid Frequencies for one treatment (combining its n-plicates)
for i in range(int(len(allTreatments)/n_plicates)): 
    treatmentList = ""
    for j in range(i*n_plicates, i*n_plicates+n_plicates):
        treatmentList+= allTreatments[j] + ','
    treatmentList = treatmentList[:-1]
    treatment = treatmentNames[i]
    print(treatmentList)
    call(["python3", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(0)])
    call(["python3", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(1)])
    call(["python3", "ObservedExpectedPlot.py", treatment])
    
treatmentList = ""
for treatment in allTreatments:
    treatmentList+= treatment
    treatmentList+=','
treatmentList = treatmentList[:-1]

#Merge_for_Deseq2.py: merges Peptide Frequnecies across all treatments into one file, for subsequent Deseq2 analysis
call(["python3", "Merge_for_Deseq2.py", experimentName, str(treatments), treatmentList])

#DESEQprocessor.r: Uses Merged Counts across all treatments for DESeq2 analysis to identify significantly enriched peptides
#                  A significantly enriched peptide is defined as p < 0.05 (Bonferroni) and log2 Fold Change > 0
subprocess.call ("DESEQprocessor.r", shell=True)


#DESEQtoFASTA.py: converts the DESEQ results into a Fasta File (sequences alone) for PCA Analysis
for i in range(len(experiments)):
    fileName1 = find_files('enrich_'+experiments[i]+'_bonf_deseq2.txt')
    print(fileName1, experiments[i])
    call(["python3", "DESEQtoFASTA.py", experiments[i], fileName1])

#PCAprocessor.r: Clusters significantly enriched peptides from Deseq using Principal Component Analysis based on Amino Acid Properties
subprocess.call ("PCAprocessor.r", shell=True)

#WeblogoProcessor.py: For each cluster (from PCA), creates a Weblogo (prereq: downloaded Weblogo as of https://github.com/WebLogo/weblogo)
treatmentList = ""
for treatment in experiments:
    treatmentList+= treatment
    treatmentList+=','
treatmentList = treatmentList[:-1]
print(treatmentList)
call(["python3", "WeblogoProcessor.py", treatmentList])
