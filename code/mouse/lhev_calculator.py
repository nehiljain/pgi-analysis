# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  29 18:00 2014

@author: nehiljain
"""

import pandas as pd
import numpy as np
import math
import csv
import os
import time
import datetime

kAnalysisName = ''

dirname = os.getcwd()
# dirname = dirname + "/Dropbox/PGI (1)/WGRanalysis"

os.chdir(dirname)

print "Code is in ", os.getcwd()

# <codecell>
lhEvOutFileName = '../'+kAnalysisName+'/Data/ProcessedData/lhEv_result'+'.csv' 
wholeEpigeneticInputFname = '../'+kAnalysisName+'/Data/ProcessedData/processed_epigenetics_raw_data.csv'

lhEvList = [["Index","Site.Start.Location","Site.End.Location","Site.Mid.Location", "Snp.Location","Lh.Snp","Snp.Name" ]]

epigeneticData = pd.read_csv(wholeEpigeneticInputFname)

print "length of the Ensemble Epigenetic Site list = ", len(epigeneticData)

testEpiData = epigeneticData[0:20]

noOfGenesSkipped = 0
snpCount = 0

# <codecell>

# def calc_lhEv(gmpl, snpl, lhsnp):
#     return float(lhsnp) * (1 -(abs(gmpl-snpl) / 500000))


# prevChromosome = 0
for epiIndex, epiRow in epigeneticData.iterrows():
    # print epiRow
    # maxLhEv = - 999999
    # maxLHSnp = - 99999
    # maxLHSnpName = ''
    # maxLhSnpLocation = - 99999
    # snpCount = 0
    # print "Gene Index and CHR", geneIndex, geneRow['Chromosome.Name']
    # if geneRow['Chromosome.Name'] is not prevChromosome:
        
    # Get the relevant lh file according to gene chromosome
    try:
        lhFilename = "../Data/RawData/alex_lhgv_calculation_files/linkage/lh_"+str(epiRow['Chromosome.Name'])+".txt"
        snpLocData = pd.read_table(lhFilename,sep="\t")
    except IOError:
        # print (lhFilename)
        noOfGenesSkipped += 1
        continue


    # print "length of the SNPs list = ", len(snpLocData)
    lowerBoundLoc = epiRow['Start.Location']
    # if lowerBoundLoc < 0:
    #     lowerBoundLoc = geneRow['Mid.Location']
    lowerIndex = np.argmin(np.abs(snpLocData['loc'] - lowerBoundLoc))
    upperIndex = np.argmin(np.abs(snpLocData['loc'] - epiRow['End.Location']))
    # upperBoundLoc = geneRow['Mid.Location'] + 500000
    
    snpLocData = snpLocData[lowerIndex:upperIndex]
    # print "length of the SNPs list = ", len(snpLocData)
    for snpIndex, snpRow in snpLocData.iterrows():
        if snpRow['loc'] >= epiRow['Start.Location'] and snpRow['loc'] <= epiRow['End.Location']:
            print(epiIndex, epiRow['Start.Location'], epiRow['End.Location'], epiRow['Mid.Location'],snpRow['loc'], snpRow['lh'], snpRow['snp_name'])
            lhEvList.append([epiIndex, epiRow['Start.Location'], epiRow['End.Location'], epiRow['Mid.Location'],snpRow['loc'], snpRow['lh'], snpRow['snp_name']])
    #         snpCount += 1
    #         currentLhEv = calc_lhEv(geneRow['Mid.Location'], snpRow['loc'], snpRow['lh'])
    #         if currentLhEv > maxLhEv:
    #             maxLhEv = currentLhEv
    #             maxLHSnp = snpRow['lh']
    #             maxLHSnpName = snpRow['snp_name']
    #             maxLhSnpLocation = snpRow['loc']
    #     elif snpRow['loc'] > upperBoundLoc:
    #         break
    # if snpCount > 0:
    #     # lhEvList.append([geneRow['Ensembl.Gene.ID'], geneRow['Mid.Location'], maxLhEv, snpCount, maxLHSnp, maxLHSnpName, maxLhSnpLocation ])
        
    # else:
    #     lhEvList.append([geneRow['Ensembl.Gene.ID'], geneRow['Mid.Location'], '', snpCount, '', '', ''])
    # # print "geneID, numberofsnps", geneRow['Ensembl.Gene.ID'], snpCount, lhEvValue
    


for chrIndex in range(1:19):
    try:
        lhFilename = "../Data/RawData/alex_lhgv_calculation_files/linkage/lh_"+str(epiRow['Chromosome.Name'])+".txt"
        snpLocData = pd.read_table(lhFilename,sep="\t")
    except IOError:
        # print (lhFilename)
        noOfGenesSkipped += 1
        continue
    snpLocData.sort(['loc'], ascending = True)
    upperEpiIndex = np.argmin(np.abs(epiRow['Chromosome.Name'] - (index+1)))
    



print "Total number of Genes Skipped = " , noOfGenesSkipped
# print "Total number of Genes Processed = " , len(geneData) - noOfGenesSkipped
############################
############################
#PRODUCING OUTPUT FILE
############################
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')

outputFile1  = open(lhEvOutFileName, "wb")
writer = csv.writer(outputFile1, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
 
for row in lhEvList: 
    writer.writerow(row)
outputFile1.close