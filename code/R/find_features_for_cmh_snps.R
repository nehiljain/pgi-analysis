library(plyr)
library(dplyr)
library(DBI)
library(RSQLite)
# library(sqldf)
library(ggplot2)
library(stringr)

setwd("~/coderepo/pgi-ngs-analysis")


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


#Data Input
featureData_chr_1 <- read.table("1_Mod_Mus.gtf", sep="\t", header = F, stringsAsFactors = F)
names(featureData_chr_1) <- c("chromosome.no", "source","feature","start","end","score","strand","frame","attributes")
snpData_chr_1 <- read.table("allSamples.chr1:1-195471971.gwas", sep="\t", header = T, stringsAsFactors = F)



featureData_chr_1 <- ddply(featureData_chr_1, .(start, end, attributes), transform, 
                                gene_id = extractGeneId(attributes),
                                gene_name = extractGeneName(attributes),
                                gene_biotype = extractGeneBiotype(attributes))


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



# findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[48000,])
# findAndStoreRelatedFeatures(featureData_chr_1,snpData_chr_1[47000,])
# 
# testSnpData <- snpData_chr_1[47100:48050,]


d_ply(snpData_chr_1, "BP", function(x) {
                          findAndStoreRelatedFeatures(featureData=featureData_chr_1,snpRow=x)
})


# Working on finding different counts
# Set Directory and Load all the libraries.

setwd("~/coderepo/pgi-ngs-analysis")
# 
library(stringr)

Df <- read.csv("SNPS_IN_FEATURES.csv", stringsAsFactors = F, row.names=NULL, header = T)

Df <- tbl_df(Df)
# Df$gene_id <- as.factor(Df$gene_id)
str(Df$gene_id)
snps_in_Gene_count_Df <- Df %.%
                            group_by(gene_id, gene_name) %.% 
                            summarise(count = n())

snps_in_features_count_Df <- Df %.%
                              group_by(feature) %.% 
                              summarise(count = n())

snps_in_sources_count_Df <- Df %.%
                              group_by(source) %.% 
                              summarise(count = n())

write.csv(snps_in_sources_count_Df, file = "snps_in_source.csv", quote=F, row.names = TRUE, col.names = FALSE)

write.csv(snps_in_features_count_Df, file = "snps_in_features_count.csv", quote=F, row.names = TRUE, col.names = FALSE)

write.csv(snps_in_Gene_count_Df, file = "snps_in_Gene_count.csv", quote=F, row.names = TRUE, col.names = FALSE)


snps_in_Gene_count_plot <- ggplot(data = snps_in_Gene_count_Df, aes(x = gene_id, y = count)) +
                              geom_bar(stat="identity", fill="#756bb1", width=.7) 

snps_in_features_count_plot <- ggplot(data = snps_in_features_count_Df, aes(x = feature, y = count)) +
                                  geom_bar(stat="identity", fill="#756bb1", width=.7) +
                                  theme(axis.text.x = element_text(angle = 90, hjust = 1))

snps_in_sources_count_plot <- ggplot(data = snps_in_sources_count_Df, aes(x = source, y = count)) +
                                geom_bar(stat="identity", fill="#756bb1", width=.7) +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(file="snps_in_Gene_count_plot.pdf", plot=snps_in_Gene_count_plot, width=12, height=8)
ggsave(file="snps_in_features_count_plot.pdf", plot=snps_in_features_count_plot, width=12, height=8)
ggsave(file="snps_in_sources_count_plot.pdf", plot=snps_in_sources_count_plot, width=12, height=8)

