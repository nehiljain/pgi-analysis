rm(list=ls())

#loading libs
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)



replaceNASnpCount <- function(snp_count, max_window_snp_count){
  if (is.na(max_window_snp_count)) {
    return(snp_count)
  }
  else return(max_window_snp_count)
}



outPrefix <- "chr1" #chromosome Name


# Name of the SNP in Gene association
inputSnpGeneFileName <- "coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv"

# Name of window analysis for that chromosome.
inputWindowAnalysisFilename <- "coderepo/pgi-ngs-analysis/data/final/Window_Analysis_Chr1.csv"


# name of outfile with all the geneInfo will be created automatically
outputFileName <- paste("coderepo/pgi-ngs-analysis/data/final/",outPrefix,"_ngs_gene_info_table",".csv", sep="")



if (file.exists(outputFileName) == TRUE) {
  file.remove(outputFileName)
}


snpGeneDf <- read.csv(inputSnpGeneFileName, stringsAsFactors = F, row.names=NULL, header = T)
snpGeneDf$start <- as.numeric(snpGeneDf$start)
snpGeneDf$end <- as.numeric(snpGeneDf$end)
snpGeneDf$snp_p_value <- as.numeric(snpGeneDf$snp_p_value)
snpGeneDf$snp_p_value <- -1 * log(snpGeneDf$snp_p_value)
snpGeneDf <- tbl_df(snpGeneDf)

windowAnalysisDf <- read.csv(inputWindowAnalysisFilename, stringsAsFactors = F, row.names=NULL, header = T)
windowAnalysisDf <- tbl_df(windowAnalysisDf)


windowAnalysisGeneDf <- windowAnalysisDf %.%
  group_by(gene_id, gene_name) %.%
  summarise(max_window_snp_count = max(snp_count))



geneInfoDf <- filter(snpGeneDf, feature == "gene") %.%
  group_by(chromosome.no, gene_id, gene_name, start, end ) %.% 
  summarise(max_cmh_neg_log = max(snp_p_value),
            snp_count = n()) %.%
  mutate(gene_length = end - start,
         average_snp_interval = (gene_length / snp_count)
  )

outputGeneInfoDf <- merge(x = geneInfoDf, y = windowAnalysisGeneDf, by = "gene_id", all.x = T )

outputGeneInfoDf$gene_name.y <- NULL


outputGeneInfoDf <- ddply(outputGeneInfoDf, .(gene_id), transform, 
              max_window_snp_count = replaceNASnpCount(snp_count, max_window_snp_count))

write.csv(outputGeneInfoDf, file = outputFileName, quote=F, row.names = FALSE)

