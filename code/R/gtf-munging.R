rm(list=ls())

library(plyr)
library(dplyr)
library(stringr)
library(rattle)

# This script loads the GTF file for mouse, creates a header, munges the strings to be made useful, converts to a dataframe and writes it to RData and CSV Files
# change the path to point to the GTF file from ensemble ftp://ftp.ensembl.org/pub/release-77/gtf/mus_musculus
# previewDf gives us idea about the structure of the file
# Header is derived from http://uswest.ensembl.org/info/website/upload/gff.html

# change the path to point to the GTF file from ensemble ftp://ftp.ensembl.org/pub/release-77/gtf/mus_musculus

# WARNING: No gene name found for line 21995 to 22008 in mouse GTF file path and url above. Replaced with "NA"
gtfFilePath <- "Downloads/Mus_musculus.GRCm38.77.gtf"

# previewDf <- read.table(file = gtfFilePath,
#                         header = FALSE,
#                         comment.char = "#",
#                         nrow = 100,
#                         na.strings = "NA",
#                         fill = TRUE,
#                         sep = "\t")


headerName <- c("chromosome_name", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

colClassNames <- c("character", "factor", "factor", "integer", "integer", "character", "character", "character", "character")
gtfData <- read.table(file = gtfFilePath,
                      header = FALSE,
                      comment.char = "#",
                      nrow = 1000,
                      na.strings = "NA",
                      fill = TRUE,
                      sep = "\t",
                      col.names = headerName)

## Reading is complete.
## Next we manipulate the attribute column to create 3 additional columns for gene id, gene name and gene biotype

removeGeneNameTag <- function(row) {
    
    s <- as.character(row)
#     print(s)
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


splitAttributes <- str_split(gtfData$attribute, "; ")

formattedAttributes <- ldply(splitAttributes, function (row) { 
    id <- removeGeneIdTag(row[grep("gene_id", row)])
    if ( length(grep("gene_name", row)) == 0 ) {
        g_name <- "NA"
        # this is done for line 21995 to 22008 in mouse GTF file path and url above.
    } else {
        g_name <- removeGeneNameTag(row[grep("gene_name", row)])   
    }
    biotype <- removeGeneBiotypeTag(row[grep("gene_biotype", row)])
    print(g_name)
    df <- data.frame(gene_id = id, gene_name = g_name, gene_biotype = biotype)
})

resultGtfData <- cbind(gtfData, formattedAttributes)

str(resultGtfData)

names(resultGtfData) <- normVarNames(names(resultGtfData))

save(resultGtfData, file = "mouse_gtf.RData")
write.csv(resultGtfData, file = "mouse_gtf.csv",
          quote = FALSE, na = "NA", row.names = FALSE)



refGeneIdData <- read.csv(, file = "Downloads/mouse_gene_list-NCBIM37.67-mm9.txt",
                          header = TRUE,
                          comment.char = "#",
                          na.strings = "NA",
                          fill = TRUE)


names(refGeneIdData) <- normVarNames(names(refGeneIdData))
save(refGeneIdData, file = "mouse_gene_list_mm9.RData")
write.csv(refGeneIdData, file = "mouse_gene_list_mm9.csv",
          quote = FALSE, na = "NA", row.names = FALSE)

write.csv(refGeneIdData[, 1], file = "test_mouse_gene_list_mm9.csv",
          quote = FALSE, na = "NA", row.names = FALSE)




