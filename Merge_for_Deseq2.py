#YuEnoch 13-11-2022
#Merge_for_Deseq2.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in all Peptide Counts for each treatment and merges all into one Table (ready to be analyzed by Deseq2), with peptides of low counts filtered out
#Changes across Experiments: alter accordingly
# 1. Filters any peptide where total counts across all treatments are below 6 (or any defined Threshold)
# 2. For experiments with 6/9/12 total Treatments (adjustable)
# 3. Mutagenesis for 5 Amino Acids

import sys
import re
import mod575
import gzip
import os
import csv
from collections import defaultdict

#Opens the necessary files for Input and Output
experimentName = sys.argv[1]
numberOfTreatments = sys.argv[2] #Example: 8
num = int(numberOfTreatments)
treatments = sys.argv[3] #Example: C1,C2,E1,E2,H1,H2,N1,N2
treatments = treatments.split(',')


combinations=open('NNK5_combinations.txt','r')	#All Possible Combinations for 5 amino acids
outfile=open(experimentName+"_merged.txt",'w')
outfileCSV=open(experimentName+"_merged.csv",'w', newline='')
writer = csv.writer(outfileCSV)
header = ["AA_seq"] + treatments
writer.writerow(header)

noreads=0
peptides=dict() #Empty Peptides Dictionary
for line in combinations:	
    line=line.strip()
    if line in peptides:	
        print('non_unique')
    else:				
        peptides[line]=0	


for i in range(num):            #Creates a Peptides Dictionary for each treatment, then goes through every Peptides File and adds to that Dictionary
    input1 = open(treatments[i]+'_peptide', 'r')
    globals()['peptide%s' % i] = peptides.copy()
    for line in input1:			
        row=line.strip().split() 	#separates into peptide, peptide count, sequence distribution
        peptide=row[0] 			
        peptide_count=row[1]
        if peptide in globals()['peptide%s' % i]:	#because peptides1 has all possible peptides, peptide must be in peptides1
            globals()['peptide%s' % i][peptide]=peptide_count
        else:
            print('sth is wrong with dictionary!')        
    input1.close()
        
for peptide in peptides:
    sum = 0
    for i in range(num):
        sum += int(globals()['peptide%s' % i][peptide])
    if sum<6: 		#Threshold of 6. If total counts for that peptide is below 6, it is not added to final result
        noreads+=1
        continue
    else:
        list = [peptide]
        for i in range(num):
            list.append(str(globals()['peptide%s' % i][peptide]))
        print ('\t'.join(list), file=outfile)
        writer.writerow(list)
outfile.close()
outfileCSV.close()

