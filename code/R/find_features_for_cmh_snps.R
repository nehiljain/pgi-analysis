rm(list=ls())

library(plyr)
library(dplyr)
library(DBI)
library(RSQLite)
library(ggplot2)
library(stringr)
# 
# setwd("~/coderepo/pgi-ngs-analysis")


# Defining the functions


extractGeneName <- function(row) {
#   print(s)
  s <- as.character(row)
  gene_name <- str_locate(s, "gene_name [a-zA-Z0-9]+")
  str_trim(str_sub(s, start = gene_name[1], end = gene_name[2]))

}


extractGeneId <- function(row) {
#   print(s)
  s <- as.character(row)
  
  gene_id <- str_locate(s, "gene_id [a-zA-Z0-9]+")
#   print(gene_id)
  str_trim(str_sub(s, start = gene_id[1], end = gene_id[2]))
}

extractGeneBiotype <- function(row) {
#   print(s)
  s <- as.character(row)
  gene_biotype<- str_locate(s, "gene_biotype [a-zA-Z0-9_]+")
  str_trim(str_sub(s, start = gene_biotype[1], end = gene_biotype[2]))

}





# create a function to get the feature file name, the gwas file name, a string suffix to be put in outfile names

featureFileName <- "data/raw/1_Mod_Mus.gtf"
gwasDirPath <- "data/raw/"
gwasFileName <- "allSamples.chr1:1-195471971.gwas"
outputPrefix <- "chr1" 
# findFeatureForCmhSnps(featureFileName, gwasFileName, outputPrefix)



findFeatureForCmhSnps <- function(featureFileName, gwasFileName, outputPrefix, outputDirPath) {

  
  #Data Input
  chrFeatureData <- read.table(featureFileName, sep="\t", header = F, stringsAsFactors = F)
  names(chrFeatureData) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")
  cmhSnpData <- read.table(gwasFileName, sep="\t", header = T, stringsAsFactors = F)
  
  cat("Dimesnions of Input Data frames", "Feature Dataset", dim(chrFeatureData), "Snps Dataset", dim(cmhSnpData), "\n\n")
  
#   separate gene attributes from the attributes column
  chrFeatureData <- ddply(chrFeatureData, .(start, end, attributes), transform, 
                          gene_id = extractGeneId(attributes),
                          gene_name = extractGeneName(attributes),
                          gene_biotype = extractGeneBiotype(attributes))
  
# Test 1   
  # cmhSnpData[48000,2]
  
  
  # testData <- filter(chrFeatureData, start <= cmhSnpData[48000,2] & end >= cmhSnpData[48000,2])
  # View(testData)
  
  #  some testing to get familiar with 
  
  
  
  dbName <- paste(outputPrefix,"_SNPS_IN_FEATURES_DATA",".sqlite", sep="")
  outFilename <- paste(outputDirPath, outputPrefix, "_SNPS_IN_FEATURES_DATA",".csv", sep="") 
  
  # Clearing the old version of the filess
  cat("dbname",dbName,"outFilename",outFilename )

  if (file.exists(outFilename) == TRUE) {
    file.remove(outFilename)
  } 
  
  db <- dbConnect(SQLite(), dbname = dbName)
  
 
  findAndStoreRelatedFeatures <- function(featureData, snpRow,db ) {
    print(snpRow$SNP)
    featuresFoundDataFrame <- filter(featureData, start <= snpRow$BP & end >= snpRow$BP)
    if (dim(featuresFoundDataFrame)[[1]] > 0) {
      featuresFoundDataFrame$snp_name <- snpRow$SNP
      featuresFoundDataFrame$snp_bp <- snpRow$BP
      featuresFoundDataFrame$snp_p_value <- snpRow$P
      
      if (dbExistsTable(conn = db, name = "SNPS_IN_FEATURES_TABLE") == TRUE) {
        dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", 
                     featuresFoundDataFrame, append = T, row.names = TRUE)
      }  else {
        dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", featuresFoundDataFrame, row.names = TRUE)
      }
      
      if (file.exists(outFilename) == TRUE) {
        write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", append = T,  col.names = FALSE, )
      }  else {
        write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", col.names = TRUE)
      }
      
    }
  }


  
#   # Test 1
#   findAndStoreRelatedFeatures(chrFeatureData,cmhSnpData[48000,])
#   # Test 2 
#   findAndStoreRelatedFeatures(featureData = chrFeatureData, snpRow = cmhSnpData[47000,])
#   # Test 3
  testSnpData <- cmhSnpData[47100:48050,]
  d_ply(testSnpData, "BP", function(x) {
    findAndStoreRelatedFeatures(featureData = chrFeatureData,snpRow = x, db)
  })

## Big Analysis 
#   d_ply(cmhSnpData, "BP", function(x) {
#     findAndStoreRelatedFeatures(featureData = chrFeatureData,snpRow = x)
#   })
  
}



cat("Starting Analysis of finding Features with SNPs within \n\n")

## Real work happening here. 

gtfDirPath = "~/coderepo/pgi-ngs-analysis/data/raw/"
gwasDirPath = "~/coderepo/pgi-ngs-analysis/data/raw/"
outputDirPath <- "data/final/"

cat("Paths \n gtfDirPath = ", gtfDirPath, " gwasDirPath = ", gwasDirPath, " outputDirPath = ", outputDirPath, "\n")

# run all chromosomes function


gtfFileNamesList <- list.files(path=gtfDirPath, pattern=".gtf", ignore.case = T)
gwasFileNamesList <- list.files(path=gwasDirPath, pattern=".gwas", ignore.case = T)

cat("List of gtfFileNames\n")
print(gtfFileNamesList)

cat("List of gwasFileNames\n")
print(gwasFileNamesList)

for ( i in 1:19) {
  cat("iteration ::",i)
  featureFileName <- gtfFileNamesList[complete.cases(
    str_locate(gtfFileNamesList, 
               pattern=paste("^",i,"_Mod_Mus.gtf",  sep="")))]


  gwasFileName <- gwasFileNamesList[complete.cases(
                                          str_locate(gwasFileNamesList, 
                                                    pattern = paste("allSamples.chr",i,":", sep="")))]
  
  cat(" featureFileName : ",featureFileName, " gwasFileName : ", gwasFileName,"\n\n")
  # print(gwasFileName)
  findFeatureForCmhSnps(featureFileName = paste(gtfDirPath,featureFileName,sep=""), 
                        gwasFileName = paste(gwasDirPath,gwasFileName,sep=""), 
                        outputPrefix = paste("chr",i, sep=""),
                        outputDirPath = outputDirPath)
  
}



