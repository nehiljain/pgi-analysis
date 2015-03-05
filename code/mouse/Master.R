
#Master File

setwd("/Users/nehiljain/Dropbox/PGI (1)/WGRanalysis")


# input a name for the analysis. eg. "wgr_is_2"  
kAnalysisName = "lhgv_133"

newFilePath1 = paste("Data", sep="/")
newFilePath4 = paste("Figures", sep="/")

newFilePath2 = paste("FinalData", sep="/")
newFilePath3 = paste("ProcessedData", sep="/")
dir.create(kAnalysisName)

setwd(kAnalysisName)
dir.create(newFilePath1)
dir.create(newFilePath4)

setwd(newFilePath1)
dir.create(newFilePath2)
dir.create(newFilePath3)


# lhgv cut-off
kLhgvLowerBound = 1.67



source('~/Dropbox/PGI (1)/WGRanalysis/Code/process_whole_mouse_reference_genome_data.r')

#Then 
source('~/Dropbox/PGI (1)/WGRanalysis/Code/merge_lhsnps.r')

#then run lhgv_calculator.py 

#then

source('~/Dropbox/PGI (1)/WGRanalysis/Code/merge_lhgv_to whole_mouse_genome.R')

#then 

source('~/Dropbox/PGI (1)/WGRanalysis/Code/david_results_from_lhgv.r')

#then david processing.py

#then
source('~/Dropbox/PGI (1)/WGRanalysis/Code/david_to_wgr_calculator.R')

#then
#clear all the variables
rm(list=ls())
source('~/Dropbox/PGI (1)/WGRanalysis/Code/top_cluster_processing.R', echo=TRUE)




