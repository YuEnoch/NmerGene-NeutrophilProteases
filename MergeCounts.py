#YuEnoch 20-02-2023
#MergeCounts.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: merge AACount or nucCount files across duplicates/triplicates for one treatment
#Changes across Experiments: please alter on parameters.txt
# 1. Location of Mutagenesis
# 2. Mutagenesis for 5 Amino Acids

import sys

number = sys.argv[1]
numberOfTreatments = sys.argv[2]
num = int(numberOfTreatments)
treatments = sys.argv[3]
treatments = treatments.split(',')
AminoOrNucleic = int(sys.argv[4]) #Merging Amino Counts (0) or Nucleic Acid Counts (1)
n_mer = int(sys.argv[5])

if AminoOrNucleic == 0:
    AAList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']    
    NNK = [2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 1, 3, 3, 2, 2, 1, 1, 1]  #default NNK codon distribution
    output = open(number+"_AACount", 'w')
    AAs = [[0 for x in range(n_mer)] for y in range(21)]
    for i in range(num):            
        input = open(treatments[i]+'_AACount', 'r')
        for line in input:
            line = line.strip().split()
            amino = line[1]
            ind = AAList.index(amino)    
            for i in range(n_mer):
                AAs[ind][i] += int(line[i+2])      
        input.close()
    for i in range(21):
        data = [number, AAList[i]]
        for j in range(n_mer):
            data.append(AAs[i][j])
        data.append(NNK[i])
        a = ""
        for k in data:
            a += str(k) + "\t"
        print(a, file=output)
        
elif AminoOrNucleic == 1:
    nucList = ['A', 'C', 'G', 'T']
    NUCs = [[0 for x in range(n_mer*3)] for y in range(4)]
    output = open(number+"_nucCount", 'w')    
    for i in range(num):            
        input = open(treatments[i]+'_nucCount', 'r')
        for line in input:
            line = line.strip().split()
            nuc = line[1]
            ind = nucList.index(nuc)    
            for i in range(n_mer*3):
                NUCs[ind][i] += int(line[i+2])      
        input.close()
    for i in range(4):
        data = [number, nucList[i]]
        for j in range(n_mer*3):
            data.append(NUCs[i][j])
        a = ""
        for k in data:
            a += str(k) + "\t"
        print(a, file=output)


output.close()
    
    
