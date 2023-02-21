#YuEnoch 20-02-2023
#getAminoAcidCount.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Peptide Data and calculates the Amino Acid Counts at each poition
#Changes across Experiments: please alter on parameters.txt
# 1. Location of Mutagenesis
# 2. Mutagenesis for 5 Amino Acids

import sys
import re
import mod575
import gzip
import os
from collections import defaultdict

#Opens the necessary files for Input and Output
number = sys.argv[1]
n_mer = int(sys.argv[2])
input = open(number+"_peptide", 'r')
AAFile = open(number+"_AACount", 'w')
AA_list=open('AA.txt','r')

AA=dict()
codons=dict()
for line in AA_list:		#AA.txt -> Reference NNK Frequencies 
    line=line.strip()
    row=line.split()	
    aa=str(row[0])	
    AA[aa]=0		
    codons[aa]=row[1]	

for i in range(n_mer):
    globals()['AA%s' % i] = AA.copy()

for line in input:
    line=line.strip()
    row=line.split()	
    peptide=row[0]		
    count=int(row[1])
    for i in range(n_mer):
        if peptide[i] in globals()['AA%s' % i]:
            globals()['AA%s' % i][peptide[i]]+=count
        else:
            globals()['AA%s' % i][peptide[i]] = count

for aa in sorted(AA):	
    data=number+'\t'+aa
    for i in range(n_mer):
        data = data + '\t' + str(globals()['AA%s' % i][aa])
    data = data + '\t' + str(codons[aa])
    print(data, file=AAFile)
    #The Frequency of Amino Acid A in Position 1, 2, 3, 4, 5, and the NNK standard

AAFile.close()
