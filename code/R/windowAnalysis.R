rm(list=ls())



library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)




outPrefix <- "chr1"
outFilename <- "coderepo/pgi-ngs-analysis/data/final/Window_Analysis_Chr1.csv"

if (file.exists(outFilename) == TRUE) {
  file.remove(outFilename=)
}

# Name iof the input file 
inputFileName <- "coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv"
Df <- read.csv(inputFileName, stringsAsFactors = F, row.names=NULL, header = T)
Df$start <- as.numeric(Df$start)
Df$end <- as.numeric(Df$end)
Df$snp_p_value <- as.numeric(Df$snp_p_value)
Df$snp_p_value <- -1 * log(Df$snp_p_value)
Df <- tbl_df(Df)



geneDfList <- dlply(filter(Df, feature == "gene"),
                    .(gene_id, gene_name))

# Constants


windowLength <- 7500
windowIncrement <- 2500

for (i  in 1:length(geneDfList)) {
  geneDf <- geneDfList[[i]]
  cat("gene name = ",geneDf$gene_name[1] ,
      "i = ", i,
      "\n")
#   str(geneDf)
  windowStart <- geneDf$start[1]
  windowEnd <- geneDf$start[1] + windowLength
  numberOfRows <- length(geneDf$gene_id)
  finalGeneDf <- data.frame(
                            snp_count = numeric(length=numberOfRows),
                            max_snp_cmh_neg_log = double(length=numberOfRows),
                            min_snp_cmh_neg_log = double(length=numberOfRows),
                            mean_snp_cmh_neg_log = double(length=numberOfRows),
                            sd_snp_cmh_neg_log = double(length=numberOfRows),
                            median_snp_cmh_neg_log = double(length=numberOfRows),
                            iqr_snp_cmh_neg_log = double(length=numberOfRows),
                            gene_id = character(length=numberOfRows),
                            gene_name = character(length=numberOfRows),
                            gene_length = integer(length=numberOfRows),
                            avg_snp_per_window = integer(length=numberOfRows), 
                            window_index = integer(length=numberOfRows),
                            stringsAsFactors=FALSE
                          )
  windowCounter <- 1
  while (windowEnd <= (floor(geneDf$end[1]))) {
    filteredGeneDf <- geneDf %.%
                    filter(snp_bp >= windowStart  & snp_bp <= windowEnd)
    if (dim(filteredGeneDf)[1] > 0) {
      windowAnalysisDf <- summarize(filteredGeneDf, snp_count = n(),
                max_snp_cmh_neg_log = max(snp_p_value),
                min_snp_cmh_neg_log = min(snp_p_value),
                mean_snp_cmh_neg_log = mean(snp_p_value),
                sd_snp_cmh_neg_log = sd(snp_p_value),
                median_snp_cmh_neg_log = median(snp_p_value),
                iqr_snp_cmh_neg_log = IQR(snp_p_value))
      windowAnalysisDf$gene_id = as.character(geneDf$gene_id[1])
      windowAnalysisDf$gene_name = as.character(geneDf$gene_name[1])
      windowAnalysisDf$gene_length = geneDf$end[1] - geneDf$start[1]
      windowAnalysisDf$window_index = windowCounter
    } else {
      windowAnalysisDf <- data.frame(
        snp_count = 0,
        max_snp_cmh_neg_log = 0,
        min_snp_cmh_neg_log = 0,
        mean_snp_cmh_neg_log = 0,
        sd_snp_cmh_neg_log = 0,
        median_snp_cmh_neg_log = 0,
        iqr_snp_cmh_neg_log = 0,
        gene_id = as.character(geneDf$gene_id[1]),
        gene_name = as.character(geneDf$gene_name[1]),
        window_index = windowCounter,
        stringsAsFactors=FALSE
      )
    }
    
    finalGeneDf[windowCounter,] =  windowAnalysisDf[1,]
    
#     cat("windowStart = ", windowStart, "windowEnd = ", windowEnd, "windowCounter = ", windowCounter, "dim filter = ", dim(windowAnalysisDf), "\n")
#     print(windowAnalysisDf[1,])                                                    
    
    windowStart <- windowStart + windowIncrement
    windowEnd <- windowStart + windowLength
    windowCounter <- windowCounter + 1
  
  }
  finalGeneDf <- finalGeneDf %.%
                    arrange(desc(snp_count), desc(max_snp_cmh_neg_log)) %.%
                    filter(window_index != 0)
                    
  print(finalGeneDf)
  windowStart <- 0
  windowEnd <- 0
  
  finalGeneDf <- finalGeneDf[,c("gene_id",
                                "gene_name",
                                "snp_count",
                                "max_snp_cmh_neg_log",
                                "min_snp_cmh_neg_log",
                                "mean_snp_cmh_neg_log",
                                "sd_snp_cmh_neg_log",
                                "median_snp_cmh_neg_log",
                                "iqr_snp_cmh_neg_log",
                                "window_index")]
  
  if (file.exists(outFilename) == TRUE) {
    write.table(finalGeneDf, file = outFilename, quote=F, sep = ",", append = T,  col.names = FALSE, row.names= FALSE)
  }  else {
    write.table(finalGeneDf, file = outFilename, quote=F, sep = ",", col.names = TRUE, row.names= FALSE)
  }
  
  rm(windowStart,windowEnd,geneDf,finalGeneDf )
}




