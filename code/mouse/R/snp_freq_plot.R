rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
library(reshape2)

genome_data_old <- read.table("/home/data/old_jn_cmh_gwas/genome.gwas", sep="\t", header = F)
genome_data <- read.table("/home/data/kacper_cmh_gwas/genome.gwas", sep="\t", header = F)
genome_data <- fread("/home/data/kacper_cmh_gwas/genome.gwas", sep="\t", sep2="auto", header=F,
      stringsAsFactors=FALSE, verbose=TRUE)
names(genome_data) <- c("chr_no", "snp_pos","ref","p_values")

# input : .gwas file tab delimited as dataframe. 
# ggplot bar plot with x values chromosme in ascending order and y snp count

get_snp_freq_plot <- function(df) {
  names(df) <- c("chr_no")
  chr_snp_freq <-  as.data.frame(table(df$chr_no))
  chr_snp_freq$Var1 <-  str_replace_all(chr_snp_freq$Var1, "chr", "")
  names(chr_snp_freq) <- c("chr_no", "snp_count")
  chr_snp_freq$chr_no <- as.integer(as.character(chr_snp_freq$chr_no))
  chr_snp_freq$chr_no[21] <- 21
  chr_snp_freq$chr_no[20] <- 20
  arrange(chr_snp_freq, chr_no)
  p1 <- ggplot(chr_snp_freq, aes(x=chr_no, y=snp_count)) + 
        geom_bar(stat="identity") +
        scale_x_continuous(breaks=1:21)
  return(p1)
}



p1 <- get_snp_freq_plot(genome_data_old)
p2 <- get_snp_freq_plot(genome_data)



p1 <- p1 + scale_y_continuous(limits=c(0,max(as.data.frame(table(genome_data$V1))[,2])))
p1 <- p1 + ylab("snp _count") + xlab("chr_no") + ggtitle("Old Analysis") + 
  theme(plot.title=element_text(family="Helvetica", face="bold", size=20) + 
          theme(axis.title=element_text(family="Helvetica", face="bold", size=20)


p2 <- p2 + scale_y_continuous(limits=c(0,max(as.data.frame(table(genome_data$V1))[,2])))
p2 <- p2 + ylab("snp _count") + xlab("chr_no") + ggtitle("New Analysis") + 
      theme(plot.title=element_text(family="Helvetica", face="bold", size=20) + 
      theme(axis.title=element_text(family="Helvetica", face="bold", size=20)
a1 <- arrangeGrob(p1,p2, nrow=2)
ggsave(filename = "/home/data/nehil_cmhtest_plots/old_vs_nwq_snp_count.png",plot=a1,scale = 1)




for(i in c("X","Y")) {
  print(i)
  chr_no <- i
  filename <- paste("/home/data/nehil_cmhtest_plots/multiple_testing_",i,".png",sep="")
  chr_df <- filter(genome_data, V1 == paste("chr",chr_no,sep=""))
  ggsave(plot=get_multiple_testing_plots(chr_df), file=filename, scale=2)
}


chr_genome_data <- ddply(genome_data, "chr_no", function(df) {
  df$chr_p_adjusted <- p.adjust(df$p_values, n=length(df$p_values), method="bonferroni")
  return(df)
})


chr_genome_data$genome_p_adjusted <- p.adjust(chr_genome_data$p_values, n=length(chr_genome_data$p_values), method="bonferroni")

names(chr_genome_data) <- c("chr_no","snp_pos","ref","p_values","chr_p_adjusted","genome_p_adjusted")
write.csv(chr_genome_data, file="/home/data/nehil_multiple_testing_csv/p_adjusted_genome.gwas", row.names=F)

plot_ready_genome_data <-melt(chr_genome_data, id.vars = c("chr_no","snp_pos","ref"))


p1 <-  ggplot(plot_ready_genome_data, aes(x = value)) + geom_histogram(bin=0.05) +
  facet_grid(variable ~ .) + 
  ggtitle("genome level")


d_ply(chr_genome_data, "chr_no", function(df) {
  i <- as.character(df$chr_no[1])
  print(i)
  plot_ready_genome_data <- melt(df, id.vars = c("chr_no","snp_pos","ref"))
  
  
  p1 <- ggplot(plot_ready_genome_data, aes(x = value)) + geom_histogram(bin=0.05) +
    facet_grid(. ~ variable) + 
    ggtitle(i)  + 
    theme(plot.title=element_text(family="Helvetica", face="bold", size=20)) + 
            theme(axis.title=element_text(family="Helvetica", face="bold", size=20))
  
  ggsave(filename = paste("/home/data/nehil_cmhtest_plots/bon_p_adjusted_",i,".png",sep=""), plot = p1, scale = 1)
})


