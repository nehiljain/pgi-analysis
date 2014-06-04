setwd("~/coderepo/pgi-ngs-analysis")
featureData_chr_1 <- read.table("1_Mod_Mus.gtf", sep="\t", header = F, stringsAsFactors = F)
names(featureData_chr_1) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")
snpData_chr_1 <- read.table("allSamples.chr1:1-195471971.gwas", sep="\t", header = F, stringsAsFactors = F)
library(dplyr)

# 
# snpData_chr_1[48000,2]


# testData <- filter(featureData_chr_1, start <= snpData_chr_1[48000,2] & end >= snpData_chr_1[48000,2])
# View(testData)

#  some testing to get familiar with 

library(DBI)
library(RSQLite)
# library(sqldf)

db <- dbConnect(SQLite(), dbname="./Test.sqlite")
# dbSendQuery(conn = db, "CREATE TABLE SNPS_IN_FEATURES (
#               ID INTEGER PRIMARY KEY   AUTOINCREMENT,
#               CHR_NO TEXT NOT NULL,
#               SOURCE TEXT NOT NULL
#             )")


# dbSendQuery(conn = db, "INSERT INTO SNPS_IN_FEATURES (CHR_NO, SOURCE) VALUES ('chr1', 'something')")
# dbWriteTable(conn = db, name = "testData1", value = testData, row.names = TRUE,)


# dbRemoveTable(conn = db, "SNPS_IN_FEATURES_TABLE")

# end of any sql testing.

# creatiion of the main function

findAndStoreRelatedFeatures <- function(featureData, snpRow ) {
  print(snpRow)
  featuresFoundDataFrame <- filter(featureData, start <= snpRow$V2 & end >= snpRow$V2)
  if (dim(featuresFoundDataFrame)[[1]] > 0) {
    featuresFoundDataFrame$snp_name <- snpRow$V3
    featuresFoundDataFrame$snp_bp <- snpRow$V2
    featuresFoundDataFrame$snp_p_value <- snpRow$V4
    print(str(featuresFoundDataFrame))
    if (dbExistsTable(conn = db, name = "SNPS_IN_FEATURES_TABLE") == TRUE) {
      dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", 
                   featuresFoundDataFrame, append = T, row.names = TRUE)
    }  else {
      dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", featuresFoundDataFrame, row.names = TRUE)
    }
  }
}


findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[48000,])
findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[47000,])


