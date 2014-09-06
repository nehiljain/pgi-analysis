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

replaceNAs <- function (row) {
    logdebug("replaceNASnpCount with params", print(row))
    if (is.na(row$max_max_snp_count)) {
        if (row$total_number_of_snps_in_gene == 1) {
            row$max_max_snp_count = row$total_number_of_snps_in_gene
            row$max_mean_snp_count = row$total_number_of_snps_in_gene
	    row$max_number_of_snps_snp_count = row$total_number_of_snps_in_gene
	    row$max_max_pvalue = row$max_cmh_neg_log
	    row$max_max_mean_pvalue = row$max_cmh_neg_log
	    row$max_mean_mean_pvalue = row$max_cmh_neg_log
            row$max_mean_pvalue = row$max_cmh_neg_log
	    row$max_number_of_snps_pvalue = row$max_cmh_neg_log
	    row$max_number_of_snps_mean_pvalue = row$max_cmh_neg_log
	   
        }
        else {
            row$max_max_snp_count = -1000
        }
    }
    return (row)
}

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

getMaxMeanRows <- function(dataframe) {
    logdebug("getMaxMeanRows with dataframe = ", str(dataframe))
    maxMaxCmhRow <- arrange(dataframe, desc(max_snp_cmh_neg_log), desc(mean_snp_cmh_neg_log))[1,]
    logdebug("maxMaxCmhRow = ", print(maxMaxCmhRow))
    maxMeanCmhRow <- arrange(dataframe, desc(mean_snp_cmh_neg_log), desc(max_snp_cmh_neg_log))[1,]
    logdebug("maxMeanCmhRow = ", print(maxMeanCmhRow))
    maxSnpCountRow <- arrange(dataframe, desc(snp_count))[1,]
    logdebug("maxSnpCountRow = ", print(maxSnpCountRow))
    combinedResult <- cbind(maxMaxCmhRow$snp_count, maxMaxCmhRow$max_snp_cmh_neg_log,  maxMaxCmhRow$mean_snp_cmh_neg_log, maxMaxCmhRow$sd_snp_cmh_neg_log,
                            maxMeanCmhRow$snp_count, maxMeanCmhRow$max_snp_cmh_neg_log,  maxMeanCmhRow$mean_snp_cmh_neg_log, maxMeanCmhRow$sd_snp_cmh_neg_log,
                            maxSnpCountRow$snp_count, maxSnpCountRow$max_snp_cmh_neg_log,  maxSnpCountRow$mean_snp_cmh_neg_log, maxSnpCountRow$sd_snp_cmh_neg_log)
    return(combinedResult)
}


getOverallGeneInfo <- function(snpGeneDf, fileName, snpsInGeneDirPath) {
    
    snpGeneDf <- filter(snpGeneDf, feature == "gene")
    logdebug("snpGeneDf", str(snpGeneDf))
   snpGeneDf <- processRefData(snpGeneDf, fileName, snpsInGeneDirPath)
    snpGeneDf <- tbl_df(snpGeneDf)
 
    geneInfoDf <- filter(snpGeneDf, feature == "gene") %.%
        group_by(gene_id) %.% 
        summarise(max_cmh_neg_log = max(snp_p_value),
                  total_number_of_snps_in_gene = n()) 
    
    
    
    snpGeneDataGeneWiseList  <- dlply(snpGeneDf,.(gene_id))
    # Get all the rows genewise 
    snpInGeneGeneWiseData <- ldply(snpGeneDataGeneWiseList, function(x) { return(x[1,1:14]) })
    
    #attach number of snps and overall max cmh value
    outputGeneInfoDf <- merge(x = geneInfoDf, y = snpInGeneGeneWiseData, by = "gene_id", all.x = T )
    
    #get gene length
    outputGeneInfoDf$gene_length = outputGeneInfoDf$end - outputGeneInfoDf$start
    
    return(outputGeneInfoDf) 
    
}




getWindowMaxRows <- function(windowAnalysisDf) {
    
    
    windowAnalysisDf <- tbl_df(windowAnalysisDf)
    windowDataGeneWiseList  <- dlply(windowAnalysisDf,.(gene_id))
    
    geneInfoFromWindowData <- ldply(windowDataGeneWiseList, function(x) { x <- getMaxMeanRows(x)
                                                                          loginfo("getMaxMeanRows returned dataframe = ", print(x))
                                                                          return(x)
    })
    
    ;
    names(geneInfoFromWindowData) <- c("gene_id", 
                                       "max_max_snp_count","max_max_pvalue","max_max_mean_pvalue","max_max_sd_pvalue",
                                       "max_mean_snp_count","max_mean_pvalue","max_mean_mean_pvalue","max_mean_sd_pvalue",
                                       "max_number_of_snps_snp_count","max_number_of_snps_pvalue","max_number_of_snps_mean_pvalue","max_number_of_snps_sd_pvalue")
    
    return(geneInfoFromWindowData)
    
}



# 
# outPrefix <- "chr1" #chromosome Name
# # Name of the SNP in Gene association
# inputSnpsInFeatureFileName <- "coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv"
# # Name of window analysis for that chromosome.
# inputWindowFileName <- "coderepo/pgi-ngs-analysis/data/sliding-window-analysis2/chr1_all-windows-sorted.csv"
# # name of outfile with all the geneInfo will be created automatically
# outputFileName <- paste("coderepo/pgi-ngs-analysis/data/final/",outPrefix,"_ngs_gene_info_table",".csv", sep="")
# if (file.exists(outputFileName) == TRUE) {
#   file.remove(outputFileName)
# }
# 
# snpGeneDf <- read.csv(inputSnpsInFeatureFileName, stringsAsFactors = F, row.names=NULL, header = T, nrow=5000)
# logdebug("snpGeneDf", str(snpGeneDf))
# #     
# #     
# outputGeneInfoDf <- getOverallGeneInfo(snpGeneDf, "inputSnpsInFeatureFileName.csv","~/coderepo/pgi-ngs-analysis/data/inputSnpGeneDirPath")
# #     
# windowAnalysisDf <- read.csv(inputWindowFileName, stringsAsFactors = F, row.names=NULL, header = T, nrow=5000)
# logdebug("windowAnalysisDf", str(windowAnalysisDf))
# geneInfoFromWindowData <- getWindowMaxRows(windowAnalysisDf)
# 
# finalDf <- merge(x = outputGeneInfoDf, y = geneInfoFromWindowData, by = "gene_id", all.x = T )
# 
# finalDf <- ddply(finalDf, .(gene_id), function(row) {
#     return(replaceNAs(row))
# } )
#     
outputDirPath = 

#  Paths required# 
inputWindowDirPath = "coderepo/pgi-ngs-analysis/data/"
inputSnpGeneDirPath <- "coderepo/pgi-ngs-analysis/data/final/"
outputGeneInfoDirPath <- "coderepo/pgi-ngs-analysis/data/"

logdebug("Paths \n inputWindowDirPath = ", print(inputWindowDirPath), " inputSnpGeneDirPath = ", print(inputSnpGeneDirPath), "\n")

# run all chromosomes function

inputSnpGeneFileNames <- list.files(path=inputSnpGeneDirPath, pattern=".csv", ignore.case = T)
inputWindowFileNames <- list.files(path=inputWindowDirPath, pattern=".csv", ignore.case = T)

cat("List of inputSnpGeneFileNames\n")
print(inputSnpGeneFileNames)



for ( i in c("X","Y")) {
    cat("iteration ::",i)
    snpsInFeatureFileName <- inputSnpGeneFileNames[complete.cases(
        str_locate(inputSnpGeneFileNames,
                   pattern=paste("chr",i,"_SNPS_IN_FEATURES_DATA",  sep="")))]
    inputSnpsInFeatureFileName <- paste(inputSnpGeneDirPath, snpsInFeatureFileName, sep="")
    
    windowsFileName <- inputWindowFileNames[complete.cases(
        str_locate(inputWindowFileNames,
                   pattern=paste("chr",i,"_all-windows-sorted",  sep="")))]
    inputWindowFileName <- paste(inputWindowDirPath, windowsFileName, sep="")
    outputFileName <- paste(outputGeneInfoDirPath,"chr",i,"_gene_info.csv", sep="")
    
    cat(" inputSnpsInFeatureFileName : ",inputSnpsInFeatureFileName, " inputWindowFileName : ", inputWindowFileName, " outputFileName : ",outputFileName,"\n\n")
    
    
    
    snpGeneDf <- read.csv(inputSnpsInFeatureFileName, stringsAsFactors = F, row.names=NULL, header = T)
    logdebug("snpGeneDf", str(snpGeneDf))
#     
#     
    outputGeneInfoDf <- getOverallGeneInfo(snpGeneDf, snpsInFeatureFileName ,inputSnpGeneDirPath)
#     
    windowAnalysisDf <- read.csv(inputWindowFileName, stringsAsFactors = F, row.names=NULL, header = T)
    logdebug("windowAnalysisDf", str(windowAnalysisDf))
    geneInfoFromWindowData <- getWindowMaxRows(windowAnalysisDf)
#     
    finalDf <- merge(x = outputGeneInfoDf, y = geneInfoFromWindowData, by = "gene_id", all.x = T )
    finalDf <- ddply(finalDf, .(gene_id), function(row) {
        return(replaceNAs(row))
    } )
    
    columnsToBeSelected <- c("chromosome.no","gene_id", "gene_name", "total_number_of_snps_in_gene","gene_length", "max_cmh_neg_log","gene_biotype","source","feature" ,"max_max_snp_count"   ,"max_max_pvalue" ,"max_max_mean_pvalue" ,"max_max_sd_pvalue"   ,"max_mean_snp_count"  ,"max_mean_pvalue" ,"max_mean_mean_pvalue","max_mean_sd_pvalue"  ,"max_number_of_snps_snp_count" ,"max_number_of_snps_pvalue" ,"max_number_of_snps_mean_pvalue" ,"max_number_of_snps_sd_pvalue", "attributes" ) 
    finalDf <- finalDf[,columnsToBeSelected]   


    write.csv(finalDf, file = outputFileName, quote=F, row.names = FALSE)
#     
  
}



# 
# 
# 
# 

