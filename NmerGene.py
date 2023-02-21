#YuEnoch 20-02-2023
#NmerGene.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: Runs All the Analyses involved, takes in 5AA Mutagenesis Data to derive a merged table (ready for Deseq2), AA_frequencies, Nucleic Acid frequencies
#         Then, it runs DESeq2 to identify significantly enriched peptides and conducts Principal Component Analyses to cluster them, which are presented as Weblogos

#         This program can accommodate any n-mer mutagenesis library at any location, for any number of experments/replicates
#         Please specify each variable in parameters.txt

#Assumptions:
# 1. Number of Replicates are same for all experiments and control. There should only be one control.
#    Example: A1, A2, A3, B1, B2, B3, N1, N2, N3 (where A/B are the experiments, N is the control)

# 2. Mutagenesis occurs as one segment between two positions (edit in parameters.txt)
#    Example: positions 55-69 inclusive (15 nucleotides, for a 5-mer) (NOTE: the position difference and n-mers should match)

# 3. Matched Seeds (the conserved sequences used) in clean_fastq.py (edit manually in clean_fastq.py)
#    NOTE: please edit manually in clean_fastq.py for your own sequencing

# 4. Quality Threshold (removes reads when any two consecutive quality scores are both below 28) (edit manually in clean_fastq.py)
# 5. Only peptides with a total count over 6 across all treatments would be analyzed in Deseq2 (edit in parameters.txt)

#Requirements:
# 1. All FastQ/Gzip files are stored in the same folder
# 2. Single-End Reads conducted in Parallel (each treatment has 2 FastQ Files, separate and not complementary)
#    Not Paired-End Reads

#NOTE: Please edit the parameters.txt file (keep the formatting same)
#      R files require packages to be installed (details in R files)
#      Ensure that the RScript Path is setup properly on your device
#           Troubleshooting: If Path is not working, replace Rscript with:
#           Windows: C:/Program Files/R/R-4.1.2/bin/Rscript.exe
#           Linux: rwrapper.executable=/usr/lib/R/bin/Rscript
#           Mac: /usr/bin/Rscript

#      Weblogo requires downloading Weblogo from https://github.com/WebLogo/weblogo

#      To run Individual Steps Separately, use # in front of lines to skip other steps


import sys
import os
import subprocess
from subprocess import call
print("starting pipeline...")

def getDetails(str, type):
    str = str.strip().split(' = ')
    str = str[1]
    if type == "s":
        return str
    elif type == "i":
        return int(str.replace(" ", ""))
    elif type == "l":
        str = str.split(', ')
        return str

parameterFile = "parameters.txt"
input=open(parameterFile,'r')

input.readline()
experimentName = getDetails(input.readline(), 's')
experiments = getDetails(input.readline(), 'l')
control = getDetails(input.readline(), 'l')

input.readline()
treatments = getDetails(input.readline(), 'i')
n_plicates = getDetails(input.readline(), 'i')  #duplicates, triplicates, etc.
treatmentNames = getDetails(input.readline(), 'l')
allTreatments = getDetails(input.readline(), 'l')

input.readline()
input.readline()
firstPosition = getDetails(input.readline(), 'i') - 1 #-1 because Python counts 0 as the 1st Position
lastPosition = getDetails(input.readline(), 'i')
n_mer = getDetails(input.readline(), 'i')
if ((lastPosition - firstPosition)/3 != n_mer):
    sys.exit("Error: Mismatch on Positions and N-Mers")
elif (len(treatmentNames) * n_plicates != treatments or len(experiments) + len(control) != len(treatmentNames) or treatments != len(allTreatments)):
    sys.exit("Error: Mismatch on Number of Treatments")
else:
    firstPosition = str(firstPosition)
    lastPosition = str(lastPosition)
    n_mer = str(n_mer)
minCount = str(getDetails(input.readline(), 'i'))

input.readline()
input.readline()
gzipChoice = getDetails(input.readline(), 's') #Are the FastQ files in Gzip? 0 = No, 1 = Yes
if gzipChoice.lower() == 'yes':
    gzipChoice = str(1)
else:
    gzipChoice = str(0)
singleParallel = getDetails(input.readline(), 's').lower()

input.readline()
input.readline()
PCAcomponents = str(getDetails(input.readline(), 'i'))
PCAclusters = str(getDetails(input.readline(), 'i'))


def fileNameFormat(num, parallel):     #FastQ File Name Format
    return allTreatments[num] + "_S" + str(num+1) + "_L00"+str(parallel)+"_R1_001.fastq"

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
    if singleParallel == 'single':
        fileName2 = 'SingleReadOnly'
    else:
        fileName2 = find_files(fileNameFormat(i, 2))
    print(name, fileName1, fileName2)
    
    call(["python3", "clean_fastq.py", name, fileName1, fileName2, gzipChoice, firstPosition, lastPosition])
    print(name, "clean_fastq.py done")
    call(["python3", "getPeptides.py", name, firstPosition, lastPosition])
    print(name, "getPeptides.py done")

    call(["python3", "getAminoAcidCount.py", name, n_mer])
    call(["python3", "getNucleicAcidCount.py", name, firstPosition, lastPosition])
    call(["python3", "ObservedExpectedPlot.py", name, n_mer])

call(["python3", "FrequencyAnalysis.py", ','.join(treatmentNames), ','.join(allTreatments)])

#Generates Observed/Expected, Nucleotide, Amino Acid Frequencies for one treatment (combining its n-plicates)
for i in range(int(len(allTreatments)/n_plicates)): 
    treatmentList = ""
    for j in range(i*n_plicates, i*n_plicates+n_plicates):
        treatmentList+= allTreatments[j] + ','
    treatmentList = treatmentList[:-1]
    treatment = treatmentNames[i]
    print(treatmentList)
    call(["python3", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(0), n_mer])
    call(["python3", "MergeCounts.py", treatment, str(n_plicates), treatmentList, str(1), n_mer])
    call(["python3", "ObservedExpectedPlot.py", treatment, n_mer])


#Merge_for_Deseq2.py: merges Peptide Frequnecies across all treatments into one file, for subsequent Deseq2 analysis
treatmentList = ','.join(allTreatments)
call(["python3", "Merge_for_Deseq2.py", experimentName, str(treatments), treatmentList, minCount])


#Variables for R Scripts
colData = ''
for treatment in experiments + control:
    for i in range(n_plicates):
        colData += treatment + ','
colData = colData[:-1]
replicates = treatmentList
experimentData = ','.join(experiments)


#DESEQprocessor.r: Uses Merged Counts across all treatments for DESeq2 analysis to identify significantly enriched peptides
#                  A significantly enriched peptide is defined as p < 0.05 (Bonferroni) and log2 Fold Change > 0
print("Conducting DESeq2 analysis:")
subprocess.call(["Rscript", os.getcwd() + "/DESEQprocessor.R", os.getcwd(), colData, replicates, experimentData, control[0], experimentName], shell=True)


#DESEQtoFASTA.py: converts the DESEQ results into a Fasta File (sequences alone) for PCA Analysis
for i in range(len(experiments)):
    fileName1 = find_files('enrich_'+experiments[i]+'_bonf_deseq2.txt')
    call(["python3", "DESEQtoFASTA.py", experiments[i], fileName1])


#PCAprocessor.r: Clusters significantly enriched peptides from Deseq using Principal Component Analysis based on Amino Acid Properties
print("Conducting PCA analysis:")
subprocess.call(["Rscript", os.getcwd() + "/PCAprocessor.r", os.getcwd(), experimentData, PCAcomponents, PCAclusters], shell=True)


#WeblogoProcessor.py: For each cluster (from PCA), creates a Weblogo
treatmentList = ""
for treatment in experiments:
    treatmentList+= treatment
    treatmentList+=','
treatmentList = treatmentList[:-1]
print(treatmentList)

print("Making Weblogos...")
call(["python3", "WeblogoProcessor.py", treatmentList])

print("Pipeline Complete")
