#YuEnoch 13-11-2022
#WeblogoProcessor.py 
#Uses Weblogo Package from https://github.com/WebLogo/weblogo

#Purpose: convert Clustering results from PCA into Weblogos 
#         also records how many peptides are in each cluster (clusterCounts.txt)

#PREREQUISITES: Weblogo Package
#To install, follow instructions on https://github.com/WebLogo/weblogo
# 1. Go to Command Prompt
# 2. Enter: pip install weblogo

import sys
from weblogo import *
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import io
import binascii
import csv    

def find(name):     #finds the file if in folder
    path = os.getcwd()
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
    return -1

def getRoot(name):     #finds the file if in folder
    path = os.getcwd()
    for root, dirs, files in os.walk(path):
        if name in files:
            return root
    return -1

outfile = open("clusterCounts.txt", 'a')

treatmentList = sys.argv[1]
treatments = treatmentList.split(',')
for i in range(len(treatments)):
    name = treatments[i]
    check = True
    cluster = 1
    
    #converts the csv file into a fasta file for Weblogo    
    while check:
        outfileName = name+'_cluster_' + str(cluster) + ".csv"
        fileLocation = find(outfileName)
        if fileLocation == -1:
            cluster = cluster - 1            
            check = False
        else:
            input = open(fileLocation, 'r')
            output1 = open(fileLocation[:-4] + ".fasta", 'w')
            output2 = open(fileLocation[:-4] + ".txt", 'w')
            
            input.readline()
            input.readline()
            count = 0
            for line in input:
                pep = line.strip().split(',')[0]
                print(">"+pep, file = output1)
                print(pep, file = output1)
                print(pep, file = output2)                
                count+=1
            print(name + ' ' + str(cluster) + ' ' + str(count), file = outfile)
            
            output1.close()
            output2.close()
            input.close()
            cluster+=1
    
    #converts the fasta file into the WebLogo
    for j in range(1,cluster+1):
        outfileName = name+'_cluster_' + str(j) + ".fasta"
        fin = open(find(outfileName), 'r')
        seqs = read_seq_data(fin)
        logodata = LogoData.from_seqs(seqs)
        logooptions = LogoOptions()
        logooptions.show_errorbars = True                
        logoformat = LogoFormat(logodata, logooptions)
        pngs = png_formatter(logodata, logoformat)
        stream = io.BytesIO(pngs)
        img = Image.open(stream)
        draw = ImageDraw.Draw(img)
        img.save(getRoot(outfileName) + '/'+name + "_Cluster_" + str(j) + "_logo.png")  
        fin.close()
        
        
outfile.close()
