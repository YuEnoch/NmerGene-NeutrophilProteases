#YuEnoch 20-02-2023
#ObservedExpectedPlot.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in Amino Acid Count Data and calculates the Observed over Expected Frequency (expected, as based on NNK)
#         to create the Observed/Expected Frequency Plot of Amino Acids
#         Additional Plots:
#         Frequency of Nucleic Acids at each Position
#         Frequency of Positions for each Amino Acid
#Changes across Experiments: please alter on parameters.txt
# 1. Location of Mutagenesis
# 2. Mutagenesis for 5 Amino Acids

import sys

NNK = [2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 1, 3, 3, 2, 2, 1, 1, 1]   #NNK Distribution, Sum is 32
#Calculates Expected Frequency of Amino Acids
ExpectedFrequency = []
for i in range(21):
    frequency = NNK[i]/32
    ExpectedFrequency.append(frequency)

number = sys.argv[1]
n_mer = int(sys.argv[2])

AAList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']
input = open(number+"_AACount", 'r')

output = open(number+"_ObsExp", 'w')
AAs = [[0 for x in range(n_mer)] for y in range(21)]
for line in input:
    line = line.strip().split()
    amino = line[1]
    ind = AAList.index(amino)    
    for i in range(n_mer):
        AAs[ind][i] = int(line[i+2])

sums = []
for i in range(21):
    a = 0
    for j in range(n_mer):
        a+=AAs[i][j]
    sums.append(a)
total = 0
for i in range(21):
    total+=sums[i]
    
#Calculates Observed Frequency of Amino Acids
ObservedFrequency = []
for i in range(21):
    frequency = sums[i]/total
    ObservedFrequency.append(frequency)
    
#Calculates Observed over Expected Frequency of Amino Acids
ObsExp = []
for i in range(21):
    frequency = ObservedFrequency[i]/ExpectedFrequency[i]
    ObsExp.append(frequency)
for i in range(21):
    print(number, AAList[i], ObsExp[i], file = output)

output2 = open(number+"_AAFreq", "w")
#Calculates Frequency of Positions for each Amino Acid
for i in range(21):
    a = number + "\t" + AAList[i] + "\t"
    for j in range(n_mer):
        frequency = AAs[i][j]/sums[i]
        a+= str(frequency) + "\t"
    print(a, file = output2)
    
#Calculates Frequency of Nucleic Acids at each Position
input = open(number+"_nucCount", 'r')
nucList = ['A', 'C', 'G', 'T']
NUC = [[0 for x in range(n_mer*3)] for y in range(4)]
for line in input:
    line = line.strip().split()
    nuc = line[1]
    ind = nucList.index(nuc)    
    for i in range(n_mer*3):
        NUC[ind][i] += int(line[i+2])      
input.close()

output2 = open(number+"_NucFreq", "w")
sums = []
for i in range(n_mer*3):
    a = 0
    for j in range(4):
        a+=NUC[j][i]
    sums.append(a)
for i in range(n_mer*3):
    a = number + "\t" + str(i) + "\t"
    for j in range(4):
        frequency = NUC[j][i]/sums[i]
        a+= str(frequency) + "\t"
    print(a, file = output2)
