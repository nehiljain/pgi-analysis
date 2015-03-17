rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)



get_multiple_testing_plots <- function(cmh_df) {
  names(cmh_df) <- c("chr_no", "snp_pos","ref","p_values")
  cmh_df <- na.omit(cmh_df)
  cmh_df$fdr_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="fdr")
  cmh_df$bon_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="bonferroni")
  # cmh_df$hommel_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="hommel")
  cmh_df$by_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="BY")
  cmh_df$bh_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="BH")
  cmh_df$hochberg_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="hochberg")
  
  cmh_df$holm_p_values <- p.adjust(cmh_df$p_value, n=length(cmh_df$p_value), method="holm")
  
  p1 <-  ggplot(cmh_df, aes(x = p_values)) + geom_histogram(bin=0.05) +
         annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
        annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, p_values < 0.05))[1]) + 
        ggtitle(levels(cmh_df$chr_no))
  p2 <- ggplot(cmh_df, aes(x = fdr_p_values)) + geom_histogram(bin=0.05) +
    annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
    annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, fdr_p_values < 0.05))[1]) + 
    ggtitle(levels(cmh_df$chr_no))
  p3 <- ggplot(cmh_df, aes(x = bon_p_values)) + geom_histogram(bin=0.05) +
    annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
    annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, bon_p_values < 0.05))[1]) + 
    ggtitle(levels(cmh_df$chr_no))
  p4 <- ggplot(cmh_df, aes(x = by_p_values)) + geom_histogram(bin=0.05) +
    annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
    annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, by_p_values < 0.05))[1]) + 
    ggtitle(levels(cmh_df$chr_no))
  p5 <- ggplot(cmh_df, aes(x = bh_p_values)) + geom_histogram(bin=0.05) +
    annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
    annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, bh_p_values < 0.05))[1]) + 
    ggtitle(levels(cmh_df$chr_no))
  p6 <- ggplot(cmh_df, aes(x = hochberg_p_values)) + geom_histogram(bin=0.05) +
    annotate("text", x = 1.5, y = 3000000, label = dim(cmh_df)[1]) +
    annotate("text", x = 1.5, y = 200000, label = dim(filter(cmh_df, hochberg_p_values < 0.05))[1]) + 
    ggtitle(levels(cmh_df$chr_no))
  a <- arrangeGrob(p1,p2,p3,p4,p5,p6,nrow=3)
  return(a)
}
