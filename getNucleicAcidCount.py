#YuEnoch 13-11-2022
#getNucleiAcidCount.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Sense sequences and calculates Nucleic Acid Frequencies at each position
#Changes across Experiments: alter accordingly
# 1. Location of Mutagenesis (here, between positions 55 and 70), the [54:69]
# 2. Mutagenesis for 5 Amino Acids

import sys
import re
import mod575
import gzip
import os
from collections import defaultdict

#Opens the necessary files for Input and Output
number = sys.argv[1]
input = open(number+"_sense", 'r')
nucFile = open(number+"_nucCount", 'w')

NUC = [[0 for x in range(15)] for y in range(4)] 


nuc_list=['A','C','G','T']
print('done')

for line in input:
    line = line.strip()
    sense=line[54:69]
    for i in range(15):
        base = sense[i]
        ind = nuc_list.index(base)
        NUC[ind][i]+=1


for i in range(4):
    Base = nuc_list[i]
    data = number+'\t'+Base+'\t'
    for j in range(15):
        data+=str(NUC[i][j]) + '\t'
#    data=number1+'\t'+nuc+'\t'+str(NUC1[nuc])+'\t'+str(NUC2[nuc])+'\t'+str(NUC3[nuc])+'\t'+str(NUC4[nuc])+'\t'+str(NUC5[nuc])+'\t'+str(NUC6[nuc])+'\t'+str(NUC7[nuc])+'\t'+str(NUC8[nuc])+'\t'+str(NUC9[nuc])+'\t'+str(NUC10[nuc])+'\t'+str(NUC11[nuc])+'\t'+str(NUC12[nuc])+'\t'+str(NUC13[nuc])+'\t'+str(NUC14[nuc])+'\t'+str(NUC15[nuc])
    
    print(data, file=nucFile)

    #1_sense   A   123   123   123   123... (15 numbers, for its frequency in each position)
    #1_sense   C   123...
    #1_sense   G   123...
    #1_sense   T   123...

nucFile.close()

