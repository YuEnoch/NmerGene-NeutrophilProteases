#YuEnoch 20-02-2023
#DESEQtoFASTA.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: converts the results from DESeq2 into Fasta files

import sys
number = sys.argv[1]
fileName = sys.argv[2]
types = sys.argv[3]
input = open(fileName, 'r')
output = open('enrich_'+number+'_'+types+'_deseq2.fasta', 'w')
input.readline()
for line in input:
    line = line.strip().split()
    print('>'+line[0], file = output)
    print(line[0], file = output)

output.close()

