#YuEnoch 13-11-2022
#getNucleiAcidCount.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Sense sequences and calculates Nucleic Acid Frequencies at each position
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
nucleotides = (int(lastPosition) - int(firstPosition))

input = open(number+"_sense", 'r')
nucFile = open(number+"_nucCount", 'w')

NUC = [[0 for x in range(nucleotides)] for y in range(4)]

nuc_list=['A','C','G','T']
print('done')

for line in input:
    line = line.strip()
    sense=line[firstPosition:lastPosition]
    for i in range(nucleotides):
        base = sense[i]
        ind = nuc_list.index(base)
        NUC[ind][i]+=1

for i in range(4):
    Base = nuc_list[i]
    data = number+'\t'+Base+'\t'
    for j in range(nucleotides):
        data+=str(NUC[i][j]) + '\t'

    print(data, file=nucFile)

    #1_sense   A   123   123   123   123... (15 numbers, for its frequency in each position)
    #1_sense   C   123...
    #1_sense   G   123...
    #1_sense   T   123...

nucFile.close()

