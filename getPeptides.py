#YuEnoch 20-02-2022
#getPeptides.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Sense Sequences of Mutagenesis Data (from clean_fastQ.py), converts into Amino Acids and counts all Peptides on a table 
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
firstPosition = int(sys.argv[2])
lastPosition = int(sys.argv[3])

fileName1 = number+'_sense'
fileName2 = number+'_peptide'
input=open(fileName1,'r')
outfile=open(fileName2,'w')

peptides=dict()
sequences=defaultdict(defaultdict)
counter=0    

for line in input:
    counter+=1
    line=line.strip()    
    sense=line[firstPosition:lastPosition]
    peptide=mod575.translate_dna(sense)
    if peptide in peptides:
        peptides[peptide]+=1
    else:
        peptides[peptide]=1

    if peptide in sequences:
        if sense in sequences[peptide]:
            sequences[peptide][sense]+=1
        else:
            sequences[peptide][sense]=1
    else:
        sequences[peptide][sense]=1

for peptide in sorted(peptides, key=peptides.get, reverse=True):
    seqlist=[]	
    for sense in sequences[peptide]:
        seqdata=str(sense)+':'+str(sequences[peptide][sense])
        seqlist.append(seqdata)
    seqs=';'.join(seqlist)
    data=peptide+'\t'+str(peptides[peptide])+'\t'+seqs
    print (data, file=outfile)

input.close()
outfile.close()
