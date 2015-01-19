rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)


removeGeneNameTag <- function(row) {
    s <- as.character(row)
    gene_name_loc <- str_locate(s, "gene_name ")
    if (!is.na(gene_name_loc[1])) {
        return(str_trim(str_sub(s, start = gene_name_loc[2])))  
    }
    return(str_trim(s))
  
}



removeGeneIdTag <- function(row) {
  s <- as.character(row)
  gene_id_loc <- str_locate(s, "gene_id ")
#   print(gene_id_loc)
  if (!is.na(gene_id_loc[1])) {
      return(str_trim(str_sub(s, start = gene_id_loc[2])))  
  }
  return(str_trim(row))

  
}

removeGeneBiotypeTag <- function(row) {
  #   print(s)
  s <- as.character(row)
  gene_biotype_loc<- str_locate(s, "gene_biotype ")
  if (!is.na(gene_biotype_loc[1])) {
      return(str_trim(str_sub(s, start = gene_biotype_loc[2])))  
  }
  return(str_trim(s))  
}



windowAnalysisPerChromosome <- function(inputFileName, outFilename) {
    if (file.exists(outFilename) == TRUE) {
        file.remove(outFilename)
    }
    cat("going to read the file\n\n")
    refDf <- read.csv(inputFileName, stringsAsFactors = F, row.names=NULL, header = T)
    refDf <- filter(refDf, feature == "gene")
    refDf$start <- as.numeric(refDf$start)
    refDf$end <- as.numeric(refDf$end)
    refDf$snp_p_value <- as.numeric(refDf$snp_p_value)
    refDf$snp_p_value <- -1 * log(refDf$snp_p_value)
    refDf <- ddply(refDf, .(snp_name, snp_p_value, gene_id, gene_name ), transform, 
                            gene_id = removeGeneIdTag(gene_id),
                            gene_name = removeGeneNameTag(gene_name),
                            gene_biotype = removeGeneBiotypeTag(gene_biotype))
    refDf$gene_id <- as.factor(refDf$gene_id)
    refDf$gene_name <- as.factor(refDf$gene_name)
    refDf$gene_biotype <- as.factor(refDf$gene_biotype)
    refDf <- tbl_df(refDf)
    
#     cat("refDf ready :: ",str(refDf))
    geneDfList <- dlply(filter(refDf, feature == "gene"),
                        .(gene_id))
    
    cat("geneDfList :: ",length(geneDfList))
    # Constants
    windowLength <- 7500
    windowIncrement <- 2500
    
    
    for (i  in 1:length(geneDfList)[]) {
        geneDf <- geneDfList[[i]]
        geneDf$gene_name <- as.character(geneDf$gene_name)
        geneDf$gene_id <- as.character(geneDf$gene_id)
        removeGeneIdTag(geneDf$gene_id)
#             cat(" <=== geneDf ---:: ",levels(geneDf$gene_id), "<---\n\n\n")
        if (as.character(geneDf$gene_name)[1] %in% c("Gm4954",  "Sema4c",  "Tmem14a")) {
            cat("gene name = ",str(as.character(geneDf$gene_name[1])),
                "i = ", i, "\n", str(geneDf))
            
        }
#         cat("gene name = ",str(as.character(geneDf$gene_name)[1]),"i = ", i, "\n")
        lastWindow <- F
        windowStart <- geneDf$start[1]
        if((floor(geneDf$end[1])) - windowStart < 7500) {
            windowLength <- (floor(geneDf$end[1])) - windowStart
        }
        windowEnd <- windowStart + windowLength
#         if (!is.na(intersect(geneDf$gene_name, c("Gm4954",  "Sema4c",  "Tmem14a")))) {
#             cat("Before:: windowStart = ", windowStart, "windowEnd = ", windowEnd, "windowlength", windowLength, "\n")
#         }
        numberOfRows <- length(geneDf$gene_id)
        finalGeneDf <- data.frame(
            snp_count = numeric(length=numberOfRows),
            max_snp_cmh_neg_log = double(length=numberOfRows),
            min_snp_cmh_neg_log = double(length=numberOfRows),
            mean_snp_cmh_neg_log = double(length=numberOfRows),
            sd_snp_cmh_neg_log = double(length=numberOfRows),
            median_snp_cmh_neg_log = double(length=numberOfRows),
            gene_id = character(length=numberOfRows),
            gene_name = character(length=numberOfRows),
            gene_length = integer(length=numberOfRows),
            window_index = integer(length=numberOfRows),
            stringsAsFactors=FALSE
        )
        windowCounter <- 1
        while ((floor(geneDf$end[1])) - windowStart > 0 & windowEnd > windowStart & !lastWindow & windowLength > 0)   {
            filteredGeneDf <- geneDf %.%
                filter(snp_bp >= windowStart  & snp_bp <= windowEnd)
            if (dim(filteredGeneDf)[1] > 0) {
                windowAnalysisDf <- summarize(filteredGeneDf, 
                                              snp_count = length(filteredGeneDf$gene_id),
                                              max_snp_cmh_neg_log = max(snp_p_value),
                                              min_snp_cmh_neg_log = min(snp_p_value),
                                              mean_snp_cmh_neg_log = mean(snp_p_value),
                                              sd_snp_cmh_neg_log = sd(snp_p_value),
                                              median_snp_cmh_neg_log = median(snp_p_value))
                windowAnalysisDf$gene_id = as.character(geneDf$gene_id[1])
                windowAnalysisDf$gene_name = as.character(geneDf$gene_name[1])
                windowAnalysisDf$gene_length = geneDf$end[1] - geneDf$start[1]
                windowAnalysisDf$window_index = windowCounter
            } else {
#                 if  (as.character(geneDf$gene_name)[1] %in% c("Gm4954",  "Sema4c",  "Tmem14a")) {
#                     cat("No SNPS gene name = ",geneDf$gene_name[1] ,
#                         "i = ", i, "\n", "NO SNPS")  
#                 }
                
                windowAnalysisDf <- data.frame(
                    snp_count = 0,
                    max_snp_cmh_neg_log = 0,
                    min_snp_cmh_neg_log = 0,
                    mean_snp_cmh_neg_log = 0,
                    sd_snp_cmh_neg_log = 0,
                    median_snp_cmh_neg_log = 0,
                    gene_id = as.character(geneDf$gene_id[1]),
                    gene_name = as.character(geneDf$gene_name[1]),
                    window_index = windowCounter,
                    stringsAsFactors=FALSE
                )
            }
            
            finalGeneDf[windowCounter,] =  windowAnalysisDf[1,]

            if  (as.character(finalGeneDf[windowCounter,"gene_name"])[1] %in% c("Gm4954",  "Sema4c",  "Tmem14a")) {
                cat("FINALLLLLL No SNPS gene name = ",str(windowAnalysisDf[1,]) ,
                    "i = ", i, "\n", "NO SNPS")  
            }
            windowStart <- windowStart + windowIncrement
            windowLength <- 7500
            if((floor(geneDf$end[1])) - windowStart < 7500) {
                windowLength <- (floor(geneDf$end[1])) - windowStart
                lastWindow <- T
            }
            windowEnd <- windowStart + windowLength
            windowCounter <- windowCounter + 1
#             cat("After:: windowStart = ", windowStart, "windowEnd = ", windowEnd, "windowCounter = ", windowCounter, "windowlength", windowLength, "\n")
            
        }

        finalGeneDf <- finalGeneDf %.%
            arrange(desc(max_snp_cmh_neg_log), desc(mean_snp_cmh_neg_log))  %.%
            filter(window_index != 0)
#         print(str(finalGeneDf))
        windowStart <- 0
        windowEnd <- 0
        windowLength <- 7500
        #   names(finalGeneDf) <- c("total_snp_count_in_gene, max_snp_cmh_neg_log, min_snp_cmh_neg_log, mean_snp_cmh_neg_log, sd_snp_cmh_neg_log, median_snp_cmh_neg_log, iqr_snp_cmh_neg_log, gene_id, gene_name,  gene_length, avg_snp_per_window, window_index")
        #   print(names(finalGeneDf))
        
        if (file.exists(outFilename) == TRUE) {
            write.table(finalGeneDf, file = outFilename, quote=F, sep = ",", append = T,  col.names = FALSE, row.names= FALSE)
        }  else {
            write.table(finalGeneDf, file = outFilename, quote=F, sep = ",", col.names = TRUE, row.names= FALSE)
        }
        
        rm(windowStart,windowEnd,geneDf,finalGeneDf )
    }
    
}



# Real Work

# required by user
# outPrefix <- "chr1"
# outFilename <- "~/coderepo/pgi-ngs-analysis/data/final/Window_Analysis_Chr1.csv"
# # Name iof the input file
# inputDirPath <- "something which holds all all files with reference features and snps found in them"
# inputFileName <- "~/coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv"
# 
# 
# # example call to function
# windowAnalysisPerChromosome(inputFileName=inputFileName, outFilename=outFilename )




cat("Starting Analysis of Window Analysis with Genes \n\n")

## Real work happening here.

outputDirPath = "/share/volatile_scratch/nehil/sliding-window-analysis/"
inputDirPath <- "/share/volatile_scratch/nehil/gene-analysis/"

cat("Paths \n outputDirPath = ", outputDirPath, " inputDirPath = ", inputDirPath, "\n")

# run all chromosomes function


inputFileNamesList <- list.files(path=inputDirPath, pattern=".csv", ignore.case = T)

cat("List of inputFileNames\n")
print(inputFileNamesList)



for ( i in 1:19) {
    cat("iteration ::",i)
    snpsInFeatureFileName <- inputFileNamesList[complete.cases(
        str_locate(inputFileNamesList,
                   pattern=paste("chr",i,"_SNPS_IN_FEATURES_DATA",  sep="")))]
    inputFileName <- paste(inputDirPath, snpsInFeatureFileName, sep="")
    outputFileName <- paste(outputDirPath,"chr",i,"_all-windows-sorted.csv", sep="")
    cat(" snpsInFeatureFileName : ",inputFileName, " output : ", outputFileName,"\n\n")
    # print(gwasFileName)
    windowAnalysisPerChromosome(inputFileName=inputFileName, outFilename=outputFileName )
}
































# 
# # Random stuff

windowDf <- read.csv("~/coderepo/pgi-ngs-analysis/data/sliding-window-analysis2/chr1_all-windows-sorted.csv")
allDf <- read.csv("~/coderepo/pgi-ngs-analysis/data/final/chr2_SNPS_IN_FEATURES_DATA.csv", row.names=NULL)


windowGeneIDList <- levels(windowDf$gene_id)
allGeneIDList <- levels(allDf$gene_id)
allGeneList <- levels(allDf$gene_name)
windowGeneList <- levels(windowDf$gene_name)

allGeneIDList <- laply(allGeneIDList, function(row) {
    print(row)
    removeGeneIdTag(row)
})
windowGeneIDList <- laply(windowGeneIDList, function(row) {
    print(row)
    removeGeneIdTag(row)
})

allGeneNameList <- laply(allGeneList, function(row) {
    print(row)
    removeGeneNameTag(row)
})
windowGeneNameList <- laply(windowGeneList, function(row) {
    print(row)
    removeGeneNameTag(row)
})

intersect(allGeneIDList, windowGeneIDList)
setdiff(allGeneIDList,windowGeneIDList)
str(windowGeneIDList)
setdiff(allGeneNameList, windowGeneList)
a <- c("a","b")
b <- c("b","c", "a")
setdiff(a,b)


someProb3 <- allDf[grep(allDf$gene_id, pattern="ENSMUSG00000025933"),]
