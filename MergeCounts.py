#YuEnoch 13-11-2022
#MergeCounts.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: merge AACount or nucCount files across duplicates/triplicates for one treatment
#Changes across Experiments: alter accordingly
# 1. Mutagenesis for 5 Amino Acids
# 2. NNK Mutagenesis

import sys

number = sys.argv[1]
numberOfTreatments = sys.argv[2] #Example: 8
num = int(numberOfTreatments)
treatments = sys.argv[3] #Example: C1,C2,E1,E2,H1,H2,N1,N2
treatments = treatments.split(',')
AminoOrNucleic = int(sys.argv[4]) #Merging Amino Counts (0) or Nucleic Acid Counts (1)

if AminoOrNucleic == 0:
    AAList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']    
    NNK = [2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 1, 3, 3, 2, 2, 1, 1, 1] 
    output = open(number+"_AACount", 'w')
    AAs = [[0 for x in range(5)] for y in range(21)] 
    for i in range(num):            
        input = open(treatments[i]+'_AACount', 'r')
        for line in input:
            line = line.strip().split()
            amino = line[1]
            ind = AAList.index(amino)    
            for i in range(5):
                AAs[ind][i] += int(line[i+2])      
        input.close()
    for i in range(21):
        data = [number, AAList[i]]
        for j in range(5):
            data.append(AAs[i][j])
        data.append(NNK[i])
        a = ""
        for k in data:
            a += str(k) + "\t"
        print(a, file=output)
        
elif AminoOrNucleic == 1:
    nucList = ['A', 'C', 'G', 'T']
    NUCs = [[0 for x in range(15)] for y in range(4)] 
    output = open(number+"_nucCount", 'w')    
    for i in range(num):            
        input = open(treatments[i]+'_nucCount', 'r')
        for line in input:
            line = line.strip().split()
            nuc = line[1]
            ind = nucList.index(nuc)    
            for i in range(15):
                NUCs[ind][i] += int(line[i+2])      
        input.close()
    for i in range(4):
        data = [number, nucList[i]]
        for j in range(15):
            data.append(NUCs[i][j])
        a = ""
        for k in data:
            a += str(k) + "\t"
        print(a, file=output)    
    
    
    
