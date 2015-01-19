#!/usr/bin/env Rscript


# find SNPS in Gene Features
# Take in Two aguments directory of
#   - SNP split chromosome wise with columns snp_id, loc
#   - GENE reference csv file
# Assumes that Gene reference file is processed and is a csv file with gene id and corresponding tables
# all paths are absolute
# example call ./

rm(list=ls())
library(plyr)
library(dplyr)
library(stringr)
library(rattle)

refGeneIdData <- NULL
lhgvDirpath <- NULL 
cmdArguments <- commandArgs(trailingOnly=TRUE)
halfWindowSize <- 500000
refGeneIdFilename <- NULL
outputDirPath <- NULL


for (i in 1:length(cmdArguments)) {
    print(paste("arg",as.character(i),"=",arguments[i]))
    outputDirPath = cmdArguments[3]
    refGeneIdFilename <- cmdArcmdArguments[2]
    lhgvDirpath <- cmdArcmdArguments[1] 
}


snpHeaderName <- c("row_name", "snp_id", "lhi", "loc")


refGeneIdData <- read.csv(file = refGeneIdFilename, header = TRUE)
names(refGeneIdData) <- normVarNames(names(refGeneIdData))
refGeneIdData$gene_mid_loc = trunc((refGeneIdData$gene_start_bp + refGeneIdData$gene_end_bp)/2)
refGeneIdData = refGeneIdData[c("ensembl_gene_id", "chromosome_name", "gene_start_bp", "gene_end_bp", "gene_mid_loc")]

chromosomeList <- c(1:19, "X", "Y")

# create a function to get the feature file name, the gwas file name, a string suffix to be put in outfile names



findFeatureFor80kSnps <- function(lhgvFileName, chromosome_number, outputDirPath) {
    
    cat("Dimesnions of Input Data frames", "Feature Dataset", dim(refGeneIdData), "Snps Dataset", dim(lhgvFileName), "\n\n")    
    outFilename <- paste(outputDirPath,"/" ,chromosome_number, "_80k_snps_in_gene_features",".csv", sep="")
    
    # Clearing the old version of the filess
    cat("outFilename",outFilename )
    
    if (file.exists(outFilename) == TRUE) {
        file.remove(outFilename)
    }
    snpColClassNames <- c("numeric","character", "numeric", "numeric")
    snpData <- read.table(file = lhgvFileName,
                              header = TRUE,
                              comment.char = "#",
                              na.strings = "NA",
                              colClasses = snpColClassNames,
                              fill = TRUE,
                              sep = "\t",
                              col.names = snpHeaderName)
    
    snpData$snp_id <- as.character(snpData$snp_id)
    snpData$lhi <- as.numeric(snpData$lhi)
    snpData$loc <- as.numeric(snpData$loc)

    findAndStoreRelatedFeatures <- function(featureData, snpRow, outFilename ) {
        print(chromosome_number)
        print(snpRow$snp_id)
        print(dim(featureData))
        featuresFoundDataFrame <- filter(featureData, 
                                         chromosome_name  ==  chromosome_number, 
                                         (gene_mid_loc - halfWindowSize) <= snpRow$loc & (gene_mid_loc + halfWindowSize) >= snpRow$loc
                                        )
        print(dim(featuresFoundDataFrame))
        if (dim(featuresFoundDataFrame)[[1]] > 0) {
            featuresFoundDataFrame$snp_id <- snpRow$snp_id
            featuresFoundDataFrame$snp_bp <- snpRow$loc
            featuresFoundDataFrame$snp_lhi <- snpRow$lhi
            if (file.exists(outFilename) == TRUE) {
                write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", append = T,  col.names = FALSE, row.names = FALSE)
            }  else {
                write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", col.names = TRUE, row.names = FALSE)
            }
            
        }
    }
    ## Big Analysis
    d_ply(snpData, .(snp_id), function(x) {
        dim(x)
        findAndStoreRelatedFeatures(featureData = filter(refGeneIdData, chromosome_name == chromosome_number), 
                                    snpRow = x, 
                                    outFilename)
    })
}
# findFeatureFor80kSnps(lhgvFileName = "~/Downloads/linkage/lh_1.txt", chromosome_number = 1, outputDirPath = "~/mel")




cat("Starting Analysis of finding Features with SNPs within \n\n")


## Real work happening here.
lhgvFilenamesList <- list.files(path=lhgvDirpath, pattern="lh_", ignore.case = T)
cat("List of lhgvFilenamesList\n")
print(lhgvFilenamesList)


# run all chromosomes function
for ( i in chromosomeList) {
    cat("iteration ::",i)
   
    lhgvFileName <- lhgvFilenamesList[complete.cases(
    str_locate(lhgvFilenamesList,
                   pattern = paste("lh_",i,".txt", sep="")))]
    lhgvFileName = paste(lhgvDirpath, lhgvFileName, sep="")
    cat(" lhgvFileName : ", lhgvFileName,"\n\n")
    # call function to start finding chromosome gene snp matches
    findFeatureFor80kSnps(lhgvFileName , chromosome_number = i, outputDirPath)
}
