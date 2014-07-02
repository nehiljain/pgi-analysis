rm(list=ls())

#loading libs
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)

outPrefix <- "chr1" #chromosome Name


# Name of the SNP in Gene association
inputFileName <- "coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv"

# Name of window analysis for that chromosome.
inputWindowAnalysisFilename <- "coderepo/pgi-ngs-analysis/data/final/Window_Analysis_Chr1.csv"




if (file.exists(outFilename) == TRUE) {
  file.remove(outFilename)
}


Df <- read.csv(inputFileName, stringsAsFactors = F, row.names=NULL, header = T)
Df$start <- as.numeric(Df$start)
Df$end <- as.numeric(Df$end)
Df$snp_p_value <- as.numeric(Df$snp_p_value)
Df$snp_p_value <- -1 * log(Df$snp_p_value)
Df <- tbl_df(Df)



metaGeneAnalysisDf <- filter(Df, feature == "gene") %.%
  group_by(chromosome.no, gene_id, gene_name, start, end ) %.% 
  summarise(max_cmh_neg_log = max(snp_p_value),
            snp_count = n()) %.%
  mutate(gene_length = end - start,
         average_snp_interval = (gene_length / snp_count)
  )

write.csv(metaGeneAnalysisDf, file = paste("data/final/",outPrefix,"gene_info_data",".csv",sep=""), quote=F, row.names = FALSE, col.names = FALSE)