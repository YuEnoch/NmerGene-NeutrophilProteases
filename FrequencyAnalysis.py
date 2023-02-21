#YuEnoch 20-02-2023
#FrequencyAnalysis.py 

#Purpose: converts the individual Frequency Analysis files into one file

import sys
import csv
treatments = sys.argv[1]
treatments = treatments.split(',')
treatmentNames = sys.argv[2] #Example: C1,C2,E1,E2,H1,H2,N1,N2
treatmentNames = treatmentNames.split(',')


AAList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']

for i in range(len(treatments)):
    input = open(treatments[i]+"_ObsExp", "r")
    globals()['list%s' % treatments[i]] = []
    for line in input:
        line = line.strip().split(' ')[2]
        globals()['list%s' % treatments[i]].append(line)
    

output = open("ObsExp_Frequencies.csv", 'w', newline = '')
writer = csv.writer(output)
header = ["Obs/Exp"] + treatments
writer.writerow(header)
for i in range(len(AAList)):
    lis = [AAList[i]]
    for j in range(len(treatments)):
        lis.append(globals()['list%s' % treatments[j]][i])
    writer.writerow(lis)

output.close()