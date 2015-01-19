rm(list=ls())

#loading libs
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(logging)



basicConfig()
setLevel('FINEST')
loginfo("Starting creation of Final gene table")

removeGeneNameTag <- function(row) {
    
    s <- as.character(row)
    print(s)
    gene_name_loc <- str_locate(s, "gene_name ")
    if (!is.na(gene_name_loc[1])) {
        return(str_trim(str_sub(s, start = gene_name_loc[2])))  
    }
    #     print(s)
    return(str_trim(s))
    
}

removeGeneIdTag <- function(row) {
    s <- as.character(row)
    #     print(s)
    gene_id_loc <- str_locate(s, "gene_id ")
    #   print(gene_id_loc)
    if (!is.na(gene_id_loc[1])) {
        return(str_trim(str_sub(s, start = gene_id_loc[2])))  
    }
    
    print(gene_id_loc)
    return(str_trim(row))
    
    
}

removeGeneBiotypeTag <- function(row) {
    
    s <- as.character(row)
    #    print(s)
    gene_biotype_loc<- str_locate(s, "gene_biotype ")
    if (!is.na(gene_biotype_loc[1])) {
        return(str_trim(str_sub(s, start = gene_biotype_loc[2])))  
    }
    return(str_trim(s))  
}

processRefData <- function(dataFrame, fileName, dirPath) {
    dataFrame <- ddply(dataFrame, .(snp_name, snp_p_value, gene_id, gene_name,start, end, attributes ), transform, 
                       gene_id = removeGeneIdTag(gene_id),
                       gene_name = removeGeneNameTag(gene_name),
                       gene_biotype = removeGeneBiotypeTag(gene_biotype))
    dataFrame$gene_id <- as.factor(dataFrame$gene_id)
    dataFrame$gene_name <- as.factor(dataFrame$gene_name)
    dataFrame$gene_biotype <- as.factor(dataFrame$gene_biotype)
    dataFrame$start <- as.numeric(dataFrame$start)
    dataFrame$end <- as.numeric(dataFrame$end)
    dataFrame$snp_p_value <- as.numeric(dataFrame$snp_p_value)
    dataFrame$snp_p_value <- -1 * log(dataFrame$snp_p_value)
    write.csv(dataFrame, file=paste(dirPath,"processed_",fileName, sep=""), quote=F)
    return(dataFrame)
}




processAndPlotGeneData <- function(snpGeneDf, fileName, snpsInGeneDirPath) {
    
    snpGeneDf <- filter(snpGeneDf, feature == "gene")
    logdebug("snpGeneDf", str(snpGeneDf))
    snpGeneDf <- processRefData(snpGeneDf, fileName, snpsInGeneDirPath)
    snpGeneDf <- tbl_df(snpGeneDf)
    d_ply(snpGeneDf, .(gene_name), function(x) {
            print(x$gene_name[1])
            plot <- ggplot(x, aes(x=snp_bp, y=snp_p_value)) +
                geom_point()
            ggsave(paste("~/pgi-ngs-analysis/figure/",x$gene_name[1],
                         "-snp-scatter-plot",fileName,".pdf",sep=""))
    })
    
    
}





#  Paths required# 
inputWindowDirPath = "~/coderepo/pgi-ngs-analysis/data/"
inputSnpGeneDirPath <- "~/coderepo/pgi-ngs-analysis/data/"
outputPlotDirPath <- "~/coderepo/pgi-ngs-analysis/figure/"

logdebug("Paths \n  inputSnpGeneDirPath = ", print(inputSnpGeneDirPath), "\n")

# run all chromosomes function

inputSnpGeneFileNames <- list.files(path=inputSnpGeneDirPath, pattern=".csv", ignore.case = T)


cat("List of inputSnpGeneFileNames\n")
print(inputSnpGeneFileNames)



for ( i in 1:length(inputSnpGeneFileNames) ) {
    cat("iteration ::",i)
    snpsInFeatureFileName <- inputSnpGeneFileNames[complete.cases(
        str_locate(inputSnpGeneFileNames,
                   pattern=paste("processed_chr4_SNPS_IN_FEATURES_DATA",  sep="")))]
    inputSnpsInFeatureFileName <- paste(inputSnpGeneDirPath, snpsInFeatureFileName, sep="")
    
    cat(" inputSnpsInFeatureFileName : ",inputSnpsInFeatureFileName, "\n\n")
    
    
    
    snpGeneDf <- read.csv(inputSnpsInFeatureFileName, stringsAsFactors = F, row.names=NULL, header = T)
    logdebug("snpGeneDf", str(snpGeneDf))
    
    
    processAndPlotGeneData(snpGeneDf, inputSnpsInFeatureFileName, inputSnpGeneDirPath)
   
}





