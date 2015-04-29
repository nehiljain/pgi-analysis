
rm(list=ls())
library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(Hmisc)

# Implementation of Heins Idea.
# The procedure below creates windows of 50 snps and calculates mean nlp for every gene.


snps_gene_df <- fread("snp_in_gene_windows.csv", sep=",",  header=T)
snps_gene_df[, nlp := -1 * log(genome_p_adjusted), by="ensembl_gene_id"]
snps_gene_df <- na.omit(g_significant)

calculate_window_mean_nlp <- function(df) {

  
  if (dim(df)[1] < 50) {
    from <- 1
    to <- dim(df)[1]
    result_df <- data.frame(mean_nlp = mean(df$nlp[from:to]),
                           chromosome_name = df$chromosome_name[1],
                           ensembl_gene_id = df$ensembl_gene_id[1],
                           window_id = "window_1")
#     cat("result df :: ",str(result_df))
  } 
  if (dim(df)[1] >= 50) {
    snp_indexes_list <- seq(50, dim(df)[1], 50)
    second_last <- tail(snp_indexes_list, n=1)
    print(snp_indexes_list)
    result_df <- ldply(snp_indexes_list, function(i) {
      from <- i-49
      to <- i
      mean_nlp <- mean(df$nlp[from:to])
      chromosome_name <- as.character(df$chromosome_name[1])
      ensembl_gene_id <- as.character(df$ensembl_gene_id[1])
      window_id <- paste("window_",i / 50)
      return(data.frame(mean_nlp, chromosome_name,ensembl_gene_id,window_id))
    })
    
#     cat("result df :: ",str(result_df))
    if (dim(df)[1] > 50) {
      from <- second_last+1
      to <- dim(df)[1]
      last_row <- data.frame(mean_nlp = mean(df$nlp[from:to]),
                             chromosome_name = as.character(df$chromosome_name[1]),
                             ensembl_gene_id = as.character(df$ensembl_gene_id[1]),
                             window_id = paste("window_",(length(snp_indexes_list) + 1)))
#       cat("Second Last",second_last,snp_indexes_list,dim(df)[1],"last row :: ",str(last_row))
      result_df <- rbind(result_df, last_row)
    }
    
  }
#   cat("result df :: ",str(result_df))
  return(result_df)
}



window_gene_df <- ddply(snps_gene_df, "ensembl_gene_id", function(df) {
#     print(str(df))
    result_df <- calculate_window_mean_nlp(df)
#     print(str(result_df))
    return(result_df)
})




