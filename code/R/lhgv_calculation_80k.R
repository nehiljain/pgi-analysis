#!/usr/bin/env Rscript

# Example Call Rscript filename.R  snpGeneDirpath  outputDirPath
# find Top SNPS in Gene Features, Calculates LHGV, Summary Stats genewise
# Take in Two aguments directory of
#   - Genes and SNps Match Table
# Assumes that  file is a csv file with gene id , multiple snps for each gene, 80k data
# all paths are absolute


rm(list=ls())
library(plyr)
library(dplyr)
library(stringr)
library(rattle)

chromosomeList <- c(1:19, "X", "Y")
snpGeneDirpath <- NULL 

cmdArguments <- commandArgs(trailingOnly=TRUE)
outputDirPath <- NULL
halfWindowSize <- 500000


for (i in 1:length(cmdArguments)) {
    print(paste("arg",as.character(i),"=",cmdArguments[i]))
    outputDirPath <- cmdArguments[2]
    snpGeneDirpath <- cmdArguments[1] 
}


# create a function to get the feature file name, the gwas file name, a string suffix to be put in outfile names






findTopSnpAndSnpStats <- function (snpGeneFilename) {
    print(snpGeneFilename)
    snpGeneDataframe = read.csv(file = snpGeneFilename, header = TRUE)
    snpGeneDataframe$snp_id <- as.character(snpGeneDataframe$snp_id)
    resultDataframe = ddply(snpGeneDataframe, .(ensembl_gene_id), function(genewiseDataframe) {

        result = as.data.frame(genewiseDataframe[1, c(1:5)])
        max_lhi_snp <- arrange(genewiseDataframe, desc(snp_lhi))[1,]
        
        result$snp_count = dim(genewiseDataframe)[1]
        result$snp_lhi_mean = mean(genewiseDataframe$snp_lhi)
        result$snp_lhi_median = median(genewiseDataframe$snp_lhi)
        result$snp_lhi_sd = sd(genewiseDataframe$snp_lhi)
        result = cbind(result, max_lhi_snp)
        return(result)
    })
    names(resultDataframe) <- c("ensembl_gene_id","chromosome_name","gene_start_bp","gene_end_bp","gene_mid_loc",
                                "genewise_snp_count","genewise_snp_lhi_mean","genewise_snp_lhi_median","genewise_snp_lhi_sd",
                                "top_snp_id","top_snp_bp","top_snp_lhi")   
    return(resultDataframe)
}


cat("Starting Analysis of finding Features with SNPs within \n\n")


## Real work happening here.
snpGeneFilenamesList <- list.files(path=snpGeneDirpath, pattern="*80k_snps_in_gene*", ignore.case = T)
cat("List of snpGeneFilenamesList\n")
print(snpGeneFilenamesList)


# run all chromosomes function
for ( i in chromosomeList) {
    print(i)
    snpGeneFilename <- snpGeneFilenamesList[complete.cases(
        str_locate(snpGeneFilenamesList,
                   pattern = paste("^",i, "_80k_",sep="")))]
    snpGeneFilename = paste(snpGeneDirpath, snpGeneFilename, sep="")
    print( snpGeneFilename)
    
    outFilename <- paste(outputDirPath,"chr_",i,"_lhgv_genes.csv", sep="")
    print(outFilename)
    
    resultDataframe = findTopSnpAndSnpStats(snpGeneFilename = snpGeneFilename)
     
    resultDataframe$lhgv = (1 - (abs(resultDataframe$gene_mid_loc - resultDataframe$top_snp_bp) / halfWindowSize)) * resultDataframe$top_snp_lhi
    
    write.csv(resultDataframe, file = outFilename)
}

