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
names(windowAnalysisDf) <- c("number_of_snps_in_window", "max_snp_cmh_neg_log", "min_snp_cmh_neg_log", "mean_snp_cmh_neg_log", "sd_snp_cmh_neg_log", "median_snp_cmh_neg_log", "iqr_snp_cmh_neg_log", "gene_id", "gene_name", "gene_length","window_index_from_start")


windowAnalysisDf <- tbl_df(windowAnalysisDf)


windowDataGeneWiseList  <- dlply(windowAnalysisDf,
                    .(gene_id, gene_name))

geneInfoFromWindowData <- ldply(windowDataGeneWiseList, function(x) {
  return(x[1,])})


geneInfoDf <- filter(snpGeneDf, feature == "gene") %.%
  group_by(chromosome.no, gene_id, gene_name, start, end ) %.% 
  summarise(max_cmh_neg_log = max(snp_p_value),
            total_number_of_snps_in_gene = n()) %.%
  mutate(gene_length = end - start,
         average_snp_location_interval_in_gene = (gene_length / total_number_of_snps_in_gene)
  )


outputGeneInfoDf <- merge(x = geneInfoDf, y = geneInfoFromWindowData, by = "gene_id", all.x = T )


outputGeneInfoDf <- ddply(outputGeneInfoDf, .(gene_id), transform, 
              max_window_snp_count = replaceNASnpCount(total_number_of_snps_in_gene, number_of_snps_in_window))
outputGeneInfoDf <- outputGeneInfoDf[ , -which(names(outputGeneInfoDf) %in% c("gene_name.y", 
                                                                              "gene_length.y",
                                                                              "avg_snp_per_window",
                                                                              ".id"))]

write.csv(outputGeneInfoDf, file = outputFileName, quote=F, row.names = FALSE)

