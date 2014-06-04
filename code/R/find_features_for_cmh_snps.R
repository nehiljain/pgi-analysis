library(dplyr)
library(DBI)
library(RSQLite)
# library(sqldf)
library(plyr)


setwd("~/coderepo/pgi-ngs-analysis")

#Data Input
featureData_chr_1 <- read.table("1_Mod_Mus.gtf", sep="\t", header = F, stringsAsFactors = F)
names(featureData_chr_1) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")
snpData_chr_1 <- read.table("allSamples.chr1:1-195471971.gwas", sep="\t", header = T, stringsAsFactors = F)




# snpData_chr_1[48000,2]


# testData <- filter(featureData_chr_1, start <= snpData_chr_1[48000,2] & end >= snpData_chr_1[48000,2])
# View(testData)

#  some testing to get familiar with 



dbName <- "./SNPS_IN_FEATURES_DATA.sqlite"
outFilename <- "./SNPS_IN_FEATURES.csv" 

# Clearing the old version of the filess
file.remove(dbName,outFilename )


db <- dbConnect(SQLite(), dbname = dbName)


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
  print(snpRow$SNP)
  featuresFoundDataFrame <- filter(featureData, start <= snpRow$BP & end >= snpRow$BP)
  if (dim(featuresFoundDataFrame)[[1]] > 0) {
    featuresFoundDataFrame$snp_name <- snpRow$SNP
    featuresFoundDataFrame$snp_bp <- snpRow$SNP
    featuresFoundDataFrame$snp_p_value <- snpRow$P
#     print(str(featuresFoundDataFrame))
    if (dbExistsTable(conn = db, name = "SNPS_IN_FEATURES_TABLE") == TRUE) {
      dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", 
                   featuresFoundDataFrame, append = T, row.names = TRUE)
    }  else {
      dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", featuresFoundDataFrame, row.names = TRUE)
    }
    if (file.exists(outFilename) == TRUE) {
      write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", append = T, row.names = TRUE, col.names = FALSE)
    }  else {
      write.table(featuresFoundDataFrame, file = outFilename, quote=F, sep = ",", row.names = TRUE, col.names = TRUE)
    }
    
  }
}


# findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[48000,])
# findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[47000,])

# testSnpData <- snpData_chr_1[48000:48050,]


d_ply(snpData_chr_1, "BP", function(x) {
                          findAndStoreRelatedFeatures(featureData=featureData_chr_1,snpRow=x)
                        })


