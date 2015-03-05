# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import numpy as np
import xlrd as xl
import math
import csv
import os
import time
import datetime

# <codecell>

dirname = os.getcwd()
# dirname = dirname + '/Dropbox/PGI (1)/WGRanalysis'
kAnalysisName = 'lhgv_167'

os.chdir(dirname)

print 'Code is in ', os.getcwd()

############################
# Importing Data
############################
headerWeightedGeneList = [['Ensembl.Gene.ID','Name.of.Included.FunctionalPathways', 
                        'Name.of.Excluded.FunctionalPathways', 'Overall.FP.Score', 'CS','CSFREQ', 'PathFreq',  
                        'PP', 'FP' ,'FuncFreq','TempNumeratorCSCalc', 'Total.Number.Of.Pathways',
                        'PPS.Of.Each.Cluster','Enrichment.Score.Each.Clutser']]
weightedGeneList = []



masterGeneList = []

gene_list = pd.read_csv('../'+kAnalysisName+'/Data/ProcessedData/top_genelist_for_wgr.csv')

print 'length of the Ensemble gene list = ', len(gene_list)

testGeneData = gene_list[0:200]

func_pathway_data = pd.read_table('../'+kAnalysisName+'/Data/david_functional_annotation_chart_report.txt', sep = '\t')

print 'length of the Ensemble gene list = ', len(func_pathway_data)

func_pathway_data = func_pathway_data.drop(['Count', '%', 'List Total', 'Pop Hits','Pop Total','Fold Enrichment','Bonferroni','Benjamini','FDR'], axis = 1)


textFilename = '../'+kAnalysisName+'/Data/david_cluster_report.txt'


cluster_data = list(csv.reader(open(textFilename, 'rb'), delimiter = '\t'))
print 'length of the cluster list = ', len(cluster_data)


############################
# Creating Basic Template
############################

for gene_index, gene_row in gene_list.iterrows():
    tempGeneRow = [ gene_row['x'], [], [], 0.0,  'NAN', 0,  0, 'NAN', 'NAN', 0, 0.0, 0, '', ''  ]
    weightedGeneList.append(tempGeneRow)
    masterGeneList.append(gene_row['x'])

testFuncData = func_pathway_data[0:10]

headerNumerOfGeneInClusterList = [['Cluster.Index','Number.of.UniqueGenes','Total.Number.of.Genes','Biological.Terms','Biological.Categories']]
numerOfGeneInClusterList = []
numberOfPathwayinClusters = 0
totalNumberOfCluster = 0


############################
# Cluster Data Processing
############################

minESValue = 1.2
maxESValue = - 99

indexOfCluster = 0

for clusterRow in cluster_data:
    # print clusterRow
    if not clusterRow:
        # print 'change of cluster \n \n'
        #do all the initialisation of for each cluster
        totalNumberOfCluster += 1
    elif 'Annotation' in clusterRow[0]:
        # print clusterRow
        clusterRow[0].strip()
        spaceIndex = clusterRow[0].rindex(' ')
        colonIndex = clusterRow[1].rindex(':')
        clusterESValue = float(clusterRow[1][colonIndex+2:])
        indexOfCluster = int(clusterRow[0][spaceIndex+1:])
        # numerOfGeneInClusterList.append([indexOfCluster, 0, 0])
        # print 'ES Score', clusterESValue, '\n'
        maxESValue = max(maxESValue, clusterESValue)
        if clusterESValue < minESValue:
            break
print minESValue, maxESValue, totalNumberOfCluster

indexOfCluster = 0
clusterGeneList = []
biologicalTermsInCluster = ''
biologialCatergoryInCluster = ''
#read only top clusters
for clusterRow in cluster_data:
    # print clusterRow
    if not clusterRow:
        print 'change of cluster ', indexOfCluster
        # do all the initialisation of for each cluster
        # print 'ES Score', indexOfCluster, clusterESValue, '\n'
        numerOfGeneInClusterList.append([indexOfCluster, len(list(set(clusterGeneList))), len(clusterGeneList), biologicalTermsInCluster, biologialCatergoryInCluster])
       
        clusterGeneList = list(set(clusterGeneList))
        
        print numerOfGeneInClusterList[indexOfCluster-1],  '\n'
        for gene in clusterGeneList:
            weightedGeneListIndex = masterGeneList.index(gene.encode('UTF-8').strip())
            weightedGeneList[weightedGeneListIndex][5] += 1
            try:
                PPScore = weightedGeneList[weightedGeneListIndex][7]/weightedGeneList[weightedGeneListIndex][6]
            except ZeroDivisionError:
                continue
                # print weightedGeneList[weightedGeneListIndex]
            weightedGeneList[weightedGeneListIndex][12] += repr(PPScore) + ','
            weightedGeneList[weightedGeneListIndex][13] += repr(clusterESValue) + ','    
            normalizedES = (clusterESValue - minESValue) / (maxESValue - minESValue)
            weightedGeneList[weightedGeneListIndex][10] += (normalizedES * PPScore)
            weightedGeneList[weightedGeneListIndex][7] = 0
            weightedGeneList[weightedGeneListIndex][6] = 0
        # perform all the calculations of a cluster for all the genes involved in the cluster before starting a new one

        clusterESValue = 0
        clusterGeneList = []
        biologicalTermsInCluster = ''
        biologialCatergoryInCluster = ''

    elif 'Annotation' in clusterRow[0]:
        # print clusterRow
        clusterRow[0].strip()
        spaceIndex = clusterRow[0].rindex(' ')
        colonIndex = clusterRow[1].rindex(':')
        clusterESValue = float(clusterRow[1][colonIndex+2:])
        indexOfCluster = int(clusterRow[0][spaceIndex+1:])
        if clusterESValue < minESValue:
            break
        print 'ES Score', indexOfCluster, clusterESValue
    elif 'ENSMUSG' in clusterRow[5]:
        print 'Pathway', clusterRow
        #each Pathway
        biologicalTermsInCluster += clusterRow[1] + ', '
        biologialCatergoryInCluster += clusterRow[0] + ', '
        numberOfPathwayinClusters += 1
        genesList = clusterRow[5].split(',')
        clusterGeneList.extend(genesList)
        pathwayPValue = float(clusterRow[4])
        logPathwayPValue = math.log(pathwayPValue) * (-1.0)
        
        for gene in genesList:
            weightedGeneListIndex = masterGeneList.index(gene.encode('UTF-8').strip())
            weightedGeneList[weightedGeneListIndex][1].extend([clusterRow[1]])
            weightedGeneList[weightedGeneListIndex][6] += 1 #Increment the freq of the gene by 1 in FuncFreq

            if weightedGeneList[weightedGeneListIndex][7] is 'NAN': # Need to check if its not yet replaced 
                # print 'Is NAN'
                weightedGeneList[weightedGeneListIndex][7] = logPathwayPValue # log PP
            
            else: 
                weightedGeneList[weightedGeneListIndex][7] += logPathwayPValue # log PP

#####################################
# Functional Pathway Processing
#####################################



countFAGene = 0
 
for funcPathIndex, funcPathRow in func_pathway_data.iterrows():
    funcPathGenelist = funcPathRow['Genes'].split(',')
    countFAGene += len(funcPathGenelist)
    # print 'funcPathIndex', funcPathIndex
    for gene in funcPathGenelist:
        weightedGeneListIndex = masterGeneList.index(gene.encode('UTF-8').strip())
        # print weightedGeneList[weightedGeneListIndex]
        if funcPathRow['Term'] not in weightedGeneList[weightedGeneListIndex][1]:
            
            weightedGeneList[weightedGeneListIndex][2].extend([funcPathRow['Term']])

            weightedGeneList[weightedGeneListIndex][9] += 1 # Increment the freq of the gene by 1 in FuncFreq
            logFP = math.log(funcPathRow['PValue']) * (-1)
            
            # weightedGeneList[weightedGeneListIndex][13] += clusterESValue + ','   @alert
            if weightedGeneList[weightedGeneListIndex][8] is 'NAN': # Need to check if its not yet replaced 
                # print 'Is NAN'
                weightedGeneList[weightedGeneListIndex][8] = logFP # log FP
                
            else: 
                weightedGeneList[weightedGeneListIndex][8] += logFP # log PP

for gene_row in weightedGeneList:
    gene_row[11] = len(gene_row[1] + gene_row[2])



############################
#PRODUCING OUTPUT FILE
############################
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')


outfileName1 = '../'+kAnalysisName+'/Data/ProcessedData/wgr_result'+'.csv'
outfileName2 = '../'+kAnalysisName+'/Data/ProcessedData/top_cluster_stats'+'.csv'

with open(outfileName1, 'wb') as f1:
    writer1 = csv.writer(f1)
    writer1.writerows(headerWeightedGeneList)
    writer1.writerows(weightedGeneList)

with open(outfileName2, 'wb') as f2:
    writer = csv.writer(f2)
    writer.writerows(headerNumerOfGeneInClusterList)
    writer.writerows(numerOfGeneInClusterList)

print 'numberOfPathwayinClusters', numberOfPathwayinClusters
print 'totalNumberOfCluster', totalNumberOfCluster
