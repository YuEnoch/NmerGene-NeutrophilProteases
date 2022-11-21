#YuEnoch 13-11-2022
#getAminoAcidCount.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Peptide Data and calculates the Amino Acid Counts at each poition
#Changes across Experiments: alter accordingly
# 1. Mutagenesis for 5 Amino Acids

import sys
import re
import mod575
import gzip
import os
from collections import defaultdict

#Opens the necessary files for Input and Output
number = sys.argv[1]
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


AA1=AA.copy()
AA2=AA.copy()
AA3=AA.copy()
AA4=AA.copy()
AA5=AA.copy()

for line in input:
    line=line.strip()
    row=line.split()	
    peptide=row[0]		
    count=int(row[1])
    if peptide[0] in AA1:
        AA1[peptide[0]]+=count	
    else:				
        AA1[peptide[0]]=count
    if peptide[1] in AA2:		
        AA2[peptide[1]]+=count
    else:
        AA2[peptide[1]]=count
    if peptide[2] in AA3:
        AA3[peptide[2]]+=count
    else:
        AA3[peptide[2]]=count
    if peptide[3] in AA4:
        AA4[peptide[3]]+=count
    else:
        AA4[peptide[3]]=count
    if peptide[4] in AA5:
        AA5[peptide[4]]+=count
    else:
        AA5[peptide[4]]=count

for aa in sorted(AA):	
    data=number+'\t'+aa+'\t'+str(AA1[aa])+'\t'+str(AA2[aa])+'\t'+str(AA3[aa])+'\t'+str(AA4[aa])+'\t'+str(AA5[aa])+'\t'+str(codons[aa])
    #Example:  99A_peptide	A	481213	476730	484370	450503	473431	4  
    #The Frequency of Amino Acid A in Position 1, 2, 3, 4, 5, and the AA.txt standard
    print(data, file=AAFile)
    
AAFile.close()