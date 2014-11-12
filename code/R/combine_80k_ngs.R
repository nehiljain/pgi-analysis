rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(xlsx)
library(rattle)
setwd("~/coderepo/pgi-ngs-analysis/data/combine_80k_ngs/")

ngsData = read.xlsx2("ngs_top_genes.xlsx", sheetName = "Sheet1")
names(ngsData) <- normVarNames(names(ngsData))

microarrayData = read.csv("wgr_annotated_top_rank_genelist.csv")
names(microarrayData) <- normVarNames(names(microarrayData))

# change name of a column 
colnames(microarrayData)[2] <- colnames(ngsData)[2]
colnames(microarrayData)[5] <- colnames(ngsData)[1]

combinedDataFrame <- merge(x = ngsData, y = microarrayData, by = "gene_id", all= FALSE)

drops <- c("x", "x_1", "chromosome_no.y")

combinedDataFrame <- combinedDataFrame[, !(names(combinedDataFrame) %in% drops )]
notFoundGeneIds <- setdiff(microarrayData$gene_id, test1$gene_id)


write.csv(combinedDataFrame, file="combined_80k_ngs_data.csv")
write.csv(notFoundGeneIds, file="genes_not_found_in_ngs.csv")



