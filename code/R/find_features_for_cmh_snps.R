library(plyr)
library(dplyr)
library(DBI)
library(RSQLite)
# library(sqldf)
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


countSnpsInSource <- function(Df, outputSuffix) {
  
  snps_in_sources_count_Df <- Df %.%
    group_by(source) %.% 
    summarise(count = n())
  write.csv(snps_in_sources_count_Df, file = paste("data/final/",outputSuffix,"snps_in_source",".csv",sep=""), quote=F, row.names = TRUE, col.names = FALSE)
  
  snps_in_sources_count_plot <- ggplot(data = snps_in_sources_count_Df, aes(x = source, y = count)) +
    geom_bar(stat="identity", fill="#756bb1", width=.7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(file=paste("figure/",outputSuffix,"snps_in_source_count",".pdf",sep=""), plot=snps_in_sources_count_plot, width=12, height=8)
}




# create a function to get the feature file name, the gwas file name, a string suffix to be put in outfile names

featureFileName <- "data/raw/1_Mod_Mus.gtf"
gwasFileName <- "data/raw/allSamples.chr1:1-195471971.gwas"
outputSuffix <- "chr1" 
# findFeatureForCmhSnps(featureFileName, gwasFileName, outputSuffix)



findFeatureForCmhSnps <- function(featureFileName, gwasFileName, outputSuffix) {

  
  #Data Input
  chrFeatureData <- read.table(featureFileName, sep="\t", header = F, stringsAsFactors = F)
  names(chrFeatureData) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")
  cmhSnpData <- read.table(gwasFileName, sep="\t", header = T, stringsAsFactors = F)
  
  
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
  
  
  
  dbName <- paste("data/final/",outputSuffix,"_SNPS_IN_FEATURES_DATA",".sqlite", sep="")
  outFilename <- paste("data/final/",outputSuffix,"_SNPS_IN_FEATURES_DATA",".csv", sep="") 
  
 
  
  
  # Clearing the old version of the filess
  
  if (file.exists(outFilename) == TRUE) {
    file.remove(outFilename)
  } 
  
  db <- dbConnect(SQLite(), dbname = dbName)
  
  
  # dbSendQuery(conn = db, "CREATE TABLE SNPS_IN_FEATURES (
  #               ID INTEGER PRIMARY KEY   AUTOINCREMENT,
  #               CHR_NO TEXT NOT NULL,
  #               SOURCE TEXT NOT NULL
  #             )")
  
  
  # dbSendQuery(conn = db, "INSERT INTO SNPS_IN_FEATURES (CHR_NO, SOURCE) VALUES ('chr1', 'something')")
  # dbWriteTable(conn = db, name = "testData1", value = testData, row.names = TRUE,)
  
  
  # dbRemoveTable(conn = db, "SNPS_IN_FEATURES_TABLE")
  
  findAndStoreRelatedFeatures <- function(featureData, snpRow ) {
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
    findAndStoreRelatedFeatures(featureData = chrFeatureData,snpRow = x)
  })

## Big Analysis 
#   d_ply(cmhSnpData, "BP", function(x) {
#     findAndStoreRelatedFeatures(featureData = chrFeatureData,snpRow = x)
#   })
  
  Df <- read.csv(outFilename, stringsAsFactors = F, row.names=NULL, header = T)
  # 
  Df <- tbl_df(Df)
  # # Df$gene_id <- as.factor(Df$gene_id)
  str(Df$gene_id)
  snps_in_Gene_count_Df <- Df %.%
    group_by(gene_id, gene_name) %.% 
    summarise(count = n())
  
  snps_in_features_count_Df <- Df %.%
    group_by(feature) %.% 
    summarise(count = n())
  # 
  

  
 
  # 
  write.csv(snps_in_features_count_Df, file = paste("data/final/",outputSuffix,"snps_in_features_count",".csv",sep=""), quote=F, row.names = TRUE, col.names = FALSE)
  # 
  write.csv(snps_in_Gene_count_Df, file = paste("data/final/",outputSuffix,"snps_in_gene_count",".csv",sep=""), quote=F, row.names = TRUE, col.names = FALSE)
  # 
  # 
  snps_in_Gene_count_plot <- ggplot(data = snps_in_Gene_count_Df, aes(x = gene_id, y = count)) +
    geom_bar(stat="identity", fill="#756bb1", width=.7) 
  # 
  snps_in_features_count_plot <- ggplot(data = snps_in_features_count_Df, aes(x = feature, y = count)) +
    geom_bar(stat="identity", fill="#756bb1", width=.7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # 
  
  
  ggsave(file=paste("figure/",outputSuffix,"snps_in_gene_count",".pdf",sep=""), plot=snps_in_Gene_count_plot, width=12, height=8)
  ggsave(file=paste("figure/",outputSuffix,"snps_in_features_count",".pdf",sep=""), plot=snps_in_features_count_plot, width=12, height=8)
 
  # 

  
}





## Real work happening here. 

gtfDirPath = "~/coderepo/pgi-ngs-analysis/data/raw/"
gwasDirPath = "~/coderepo/pgi-ngs-analysis/data/raw/"


# run all chromosomes function


gtfFileNamesList <- list.files(path=gtfDirPath, pattern=".gtf", ignore.case = T)
gwasFileNamesList <- list.files(path=gwasDirPath, pattern=".gwas", ignore.case = T)


for ( i in 1:19) {
  print(i)
  featureFileName <- gtfFileNamesList[complete.cases(
    str_locate(gtfFileNamesList, 
               pattern=paste("^",i,"_Mod_Mus.gtf",  sep="")))]


  gwasFileName <- gwasFileNamesList[complete.cases(
                                          str_locate(gwasFileNamesList, 
                                                    pattern = paste("allSamples.chr",i,":", sep="")))]
  
  print(featureFileName)
  print(gwasFileName)
  findFeatureForCmhSnps(featureFileName = paste(gtfDirPath,featureFileName,sep=""), 
                        gwasFileName = paste(gwasDirPath,gwasFileName,sep=""), 
                        outputSuffix = paste("chr",i, sep="_"))
  
}



