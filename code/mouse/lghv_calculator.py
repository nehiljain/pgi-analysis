# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 22:43:23 2014

@author: nehiljain
"""

import pandas as pd
import numpy as np
import math
import csv
import os
import time
import datetime

kAnalysisName = 'lhgv_167'

dirname = os.getcwd()
# dirname = dirname + "/Dropbox/PGI (1)/WGRanalysis"

os.chdir(dirname)

print "Code is in ", os.getcwd()

# <codecell>
lhgvOutFileName = '../'+kAnalysisName+'/Data/ProcessedData/lhgv_result'+'.csv' 
wholeGenomeforlhSNPCalcInputFname = '../'+kAnalysisName+'/Data/ProcessedData/whole_mouse_genome_for_lhSnp_calc.csv'

lhgvList = [["Ensembl.Gene.ID","Gene.Mid.Location","LHGV","No.Of.Snps","Top.Lh.Snp","Top.Snp.Name", "Top.Snp.Location"]]

geneData = pd.read_csv(wholeGenomeforlhSNPCalcInputFname)

print "length of the Ensemble gene list = ", len(geneData)

testGeneData = geneData[0:2000]

noOfGenesSkipped = 0
snpCount = 0

# <codecell>

def calc_lhgv(gmpl, snpl, lhsnp):
    return float(lhsnp) * (1 -(abs(gmpl-snpl) / 500000))


prevChromosome = 0
for geneIndex, geneRow in geneData.iterrows():
    print "."
    maxLhgv = - 999999
    maxLHSnp = - 99999
    maxLHSnpName = ''
    maxLhSnpLocation = - 99999
    snpCount = 0
    # print "Gene Index and CHR", geneIndex, geneRow['Chromosome.Name']
    # if geneRow['Chromosome.Name'] is not prevChromosome:
        
    # Get the relevant lh file according to gene chromosome
    try:
        lhFilename = "../Data/RawData/alex_lhgv_calculation_files/linkage/lh_"+geneRow['Chromosome.Name']+".txt"
        snpLocData = pd.read_table(lhFilename,sep="\t")
    except IOError:
        print ("ERROR :: Gene Index", geneIndex, "and Chromosome: ", geneRow['Chromosome.Name'])
        noOfGenesSkipped += 1
        continue


    # print "length of the SNPs list = ", len(snpLocData)
    lowerBoundLoc = geneRow['Mid.Location'] - 500000
    if lowerBoundLoc < 0:
        lowerBoundLoc = geneRow['Mid.Location']
    lowerIndex = np.argmin(np.abs(snpLocData['loc'] - lowerBoundLoc))
    upperBoundLoc = geneRow['Mid.Location'] + 500000
    
    snpLocData = snpLocData[lowerIndex:]
    # print "length of the SNPs list = ", len(snpLocData)
    for snpIndex, snpRow in snpLocData.iterrows():
        if snpRow['loc'] >= lowerBoundLoc and snpRow['loc'] <= upperBoundLoc:
            snpCount += 1
            currentLhgv = calc_lhgv(geneRow['Mid.Location'], snpRow['loc'], snpRow['lh'])
            if currentLhgv > maxLhgv:
                maxLhgv = currentLhgv
                maxLHSnp = snpRow['lh']
                maxLHSnpName = snpRow['snp_name']
                maxLhSnpLocation = snpRow['loc']
        elif snpRow['loc'] > upperBoundLoc:
            break
    if snpCount > 0:
        lhgvList.append([geneRow['Ensembl.Gene.ID'], geneRow['Mid.Location'], maxLhgv, snpCount, maxLHSnp, maxLHSnpName, maxLhSnpLocation ])
    else:
        lhgvList.append([geneRow['Ensembl.Gene.ID'], geneRow['Mid.Location'], '', snpCount, '', '', ''])
    # print "geneID, numberofsnps", geneRow['Ensembl.Gene.ID'], snpCount, lhgvValue
    prevChromosome = geneRow['Chromosome.Name']

print "Total number of Genes Skipped = " , noOfGenesSkipped
print "Total number of Genes Processed = " , len(geneData) - noOfGenesSkipped
############################
############################
#PRODUCING OUTPUT FILE
############################
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')

outputFile1  = open(lhgvOutFileName, "wb")
writer = csv.writer(outputFile1, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
 
for row in lhgvList: 
    writer.writerow(row)
outputFile1.close