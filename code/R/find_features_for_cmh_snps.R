setwd("~/coderepo/pgi-ngs-analysis")
featureData_chr_1 <- read.table("1_Mod_Mus.gtf", sep="\t", header = F, stringsAsFactors = F)
snpData_chr_1 <- read.table("allSamples.chr1:1-195471971.gwas", sep="\t", header = F, stringsAsFactors = F)
library(dplyr)
filter(featureData_chr_1, start <= snpData_chr_1[48000,2] & end <= snpData_chr_1[48000,2])
snpData_chr_1[48000,2]
names(featureData_chr_1) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")

testData <- filter(featureData_chr_1, start <= snpData_chr_1[48000,2] & end >= snpData_chr_1[48000,2])
View(testData)

#  some testing to get familiar with 

library(DBI)
library(RSQLite)
library(sqldf)

db <- dbConnect(SQLite(), dbname="./RemoveMe.sqlite")
dbSendQuery(conn = db, "CREATE TABLE SNPS_IN_FEATURES (
              ID INTEGER PRIMARY KEY   AUTOINCREMENT,
              CHR_NO TEXT NOT NULL,
              SOURCE TEXT NOT NULL
            )")


dbSendQuery(conn = db, "INSERT INTO SNPS_IN_FEATURES (CHR_NO, SOURCE) VALUES ('chr1', 'something')")
dbWriteTable(conn = db, name = "testData1", value = testData, row.names = TRUE,)


if(dbExistsTable(conn = db, name = "SNPS_IN_FEATURES_TABLE") == TRUE) {
  dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", testData, append = T, row.names = TRUE)
}  else {
    dbWriteTable(conn = db,  name = "SNPS_IN_FEATURES_TABLE", testData, row.names = TRUE)
}
