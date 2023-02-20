#YuEnoch 20-02-2023
#clean_fastQ.py 
#Based off code from Dr. Kart Tomberg -> https://github.com/tombergk/NNK_VWF73/

#Purpose: takes in raw FastQ data (single-end read) and filters sequencing reads by Quality and Matched Seeds, to obtain good Sense sequences (and Junk data)
#Changes across Experiments: please alter on parameters.txt
# 1. Location of Mutagenesis
# 2. Mutagenesis for 5 Amino Acids
# 3. Sequencing Reads are done in Parallel (two separate reads done at same time)
#    If done in Single, please alter on parameters.txt
#    This program only takes One-End Reads, not Paired-End Reads

#Manual Changes: please alter within this Python Script
# 1. Matched Seeds (the conserved sequences used)
# 2. Quality Threshold (here, removes reads when any two consecutive quality scores are both below 28)

import sys
import re
import mod575
import gzip
import os
from collections import defaultdict

def get_read_set(file1):
    '''Expects to open file as input. Returns 4 lines of data, four from
each file. Return characters are stripped from each line.
'''
    if gzipChoice == 0:
        while True:
            line1 = str(file1.readline().strip())
            line2 = str(file1.readline().strip())
            line3 = str(file1.readline().strip())
            line4 = str(file1.readline().strip())
            if (line1 and line2 and line3 and line4):
                yield line1, line2, line3, line4
            else:
                break
    elif gzipChoice == 1:
        while True:
            line1 = str(file1.readline().strip())[2:-1] 
            line2 = str(file1.readline().strip())[2:-1]
            line3 = str(file1.readline().strip())[2:-1]
            line4 = str(file1.readline().strip())[2:-1]
            if (line1 and line2 and line3 and line4):
                yield line1, line2, line3, line4
            else:
                break

#Opens the necessary files for Input and Output
number = sys.argv[1]
fileName1 = sys.argv[2]
fileName2 = sys.argv[3]
gzipChoice = int(sys.argv[4])
firstPosition = int(sys.argv[5])
lastPosition = int(sys.argv[6])

if gzipChoice == 0:
    file1 = open(fileName1, 'r')
    file2 = open(fileName2, 'r')
elif gzipChoice == 1:
    file1 = gzip.open(fileName1, 'r')
    file2 = gzip.open(fileName2, 'r')

singleRead = False
if fileName2 == "SingleReadOnly":
    singleRead = True

good1=open(number+"_good_fastq",'w') 
sense1=open(number+"_sense",'w')
strand1=open(number+"_strand",'w')
junk1=open(number+"_junk.fastq", 'w')
stats=open('stats_clean_fastq.txt','a')

#Counters to keep track of Stats and Process
counter=0                   
good_reads=0
seed1=0
seed2=0
seed3=0
seed4=0
seed_unmatched=0
bad_quality=0
sense_antisense_unmatched=0

for pairedRead in get_read_set(file1):
    counter+=1 
    row1=pairedRead[0]
    row2=pairedRead[1]
    row3=pairedRead[2]
    row4=pairedRead[3]

    #Checks the read for matching of at least one of three seeds of 8bp. If no match, goes to junk
    if re.match(r'^CAGGAGGA[ATCGN]+',row2):
        m=re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*',row2)
        sense=m.group(1)
        sense_quality=row4[firstPosition:lastPosition]
        sense_row=row2
        seed1+=1
    elif re.match(r'^[ATCGN]{12}CAACAAAT[ATCGN]+',row2):
        m=re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*',row2)
        sense=m.group(1)
        sense_quality=row4[firstPosition:lastPosition]
        sense_row=row2
        seed2+=1
    elif re.match(r'^[ATCGN]{20}GGGTGGTG[ATCGN]+',row2):
        m=re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*',row2)
        sense=m.group(1)
        sense_quality=row4[firstPosition:lastPosition]
        sense_row=row2
        seed3+=1
    else:
        seed_unmatched+=1
        print(row1, file=junk1)
        print(row2, file=junk1)
        print(row3, file=junk1)
        print(row4, file=junk1)
        continue

    #Quality Check (Scale of 0-41): removes any read where two consecutive bases have a quality below 28
    poorQuality = False
    poor_quality = 0
    for i in range(len(sense_quality)-1):
        pos1 = sense_quality[i]
        pos2 = sense_quality[i+1]
        score1=int(ord(str(pos1)))-33
        score2=int(ord(str(pos2)))-33
        if score1 < 28 and score2 < 28:
            poorQuality = True

    if poorQuality:
        poor_quality+=1
        print(row1, file=junk1)
        print(row2, file=junk1)
        print(row3, file=junk1)
        print(row4, file=junk1)
        continue

    else:           #Those who pass the Quality Filter are kept
        good_reads+=1
        print(row1, file=good1)
        print(row2, file=good1)
        print(row3, file=good1)
        print(row4, file=good1)
        print(sense_row, file=sense1)
        print(sense_row[firstPosition:lastPosition], file=strand1)

if not singleRead:
    for pairedRead in get_read_set(file2):
        counter += 1
        row1 = pairedRead[0]
        row2 = pairedRead[1]
        row3 = pairedRead[2]
        row4 = pairedRead[3]

        # Checks the read for matching of at least one of three seeds of 8bp. If no match, goes to junk
        if re.match(r'^CAGGAGGA[ATCGN]+', row2):
            m = re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*', row2)
            sense = m.group(1)
            sense_quality = row4[firstPosition:lastPosition]
            sense_row = row2
            seed1 += 1
        elif re.match(r'^[ATCGN]{12}CAACAAAT[ATCGN]+', row2):
            m = re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*', row2)
            sense = m.group(1)
            sense_quality = row4[firstPosition:lastPosition]
            sense_row = row2
            seed2 += 1
        elif re.match(r'^[ATCGN]{20}GGGTGGTG[ATCGN]+', row2):
            m = re.match(r'^[ATCGN]{50}([ATCGN]{15})[ATCGN]*', row2)
            sense = m.group(1)
            sense_quality = row4[firstPosition:lastPosition]
            sense_row = row2
            seed3 += 1
        else:
            seed_unmatched += 1
            print(row1, file=junk1)
            print(row2, file=junk1)
            print(row3, file=junk1)
            print(row4, file=junk1)
            continue

        # Quality Check (Scale of 0-41): removes any read where two consecutive bases have a quality below 28
        poorQuality = False
        poor_quality = 0
        for i in range(len(sense_quality) - 1):
            pos1 = sense_quality[i]
            pos2 = sense_quality[i + 1]
            score1 = int(ord(str(pos1))) - 33
            score2 = int(ord(str(pos2))) - 33
            if score1 < 28 and score2 < 28:
                poorQuality = True

        if poorQuality:
            poor_quality += 1
            print(row1, file=junk1)
            print(row2, file=junk1)
            print(row3, file=junk1)
            print(row4, file=junk1)
            continue

        else:  # Those who pass the Quality Filter are kept
            good_reads += 1
            print(row1, file=good1)
            print(row2, file=good1)
            print(row3, file=good1)
            print(row4, file=good1)
            print(sense_row, file=sense1)
            print(sense_row[firstPosition:lastPosition], file=strand1)

        
        
print(number,counter,good_reads,seed1,seed2,seed3,seed_unmatched,bad_quality, file=stats)

file1.close()
good1.close()
junk1.close()
sense1.close()
stats.close()
strand1.close()
