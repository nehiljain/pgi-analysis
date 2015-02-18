library(RDAVIDWebService)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

kLhgvLowerBound = 1

mergedLHGVOutputFname = "ProcessedData/merged_lhgv_results_whole_mouse_genome.csv" 

# Read Genelist 
wholeGeneList <-  read.csv(mergedLHGVOutputFname, 
                           header=TRUE, comment.char="", nrows=40000)

#Descending Order
wholeGeneList <- wholeGeneList[with(wholeGeneList, order(-wholeGeneList$LHGV)), ]


david <- DAVIDWebService$new(email = 'nehil.jain@mail.mcgill.ca')
connect(david)
idTypes <- getIdTypes(david)


missingValue <- is.na(wholeGeneList$LHGV)
wholeGeneList[missingValue,"LHGV"] <- - 99.99
# wholeGeneList <- wholeGeneList[wholeGeneList$Gene.Biotype == "Amino_Acid_Coding",]
topLhGVBool <- wholeGeneList$LHGV > kLhgvLowerBound

topGeneList <- wholeGeneList[topLhGVBool,"Ensembl.Gene.ID"]

numberOfIterations = 100
idIncrements = 32
metaDataFrame = data.frame(matrix(vector(), numberOfIterations, 7, dimnames=list(c(), c("Number.Of.Genes",
                                                                                        "Number.Of.Clusters",
                                                        
                                                               "Max.ES","Min.ES",
                                                               "Median.ES","SD.ES","Mean.ES"))), stringsAsFactors=F)
# metaDataFrame[1,] = c(11, "ruyiuyt", "hkgkhgk",9,0,9,9,9)

upperBoundIndex = 3200+idIncrements

for (i in 1:numberOfIterations){
  print(i)
  geneList <- topGeneList[1:upperBoundIndex]
  addListResult <- addList(david, geneList, idType="ENSEMBL_GENE_ID", listName = paste("mouseGeneList",i,sep=""), listType=c("Gene"))
  clusterReport <- getClusterReport(david)
  higherEnrichment <- (enrichment(clusterReport)>1.0)
  listES <- summary(clusterReport)[higherEnrichment,c("Enrichment")]
  metaDataFrame[i,c("Number.Of.Genes")] <- upperBoundIndex
  metaDataFrame[i,c("Number.Of.Clusters")] <- length(listES)
  metaDataFrame[i,c("Mean.ES")] <- mean(listES)
  metaDataFrame[i, c("Max.ES")] <- max(listES)
  metaDataFrame[i, c("Min.ES")] <- min(listES)
  metaDataFrame[i, c("SD.ES")] <- sd(listES)
  metaDataFrame[i, c("Median.ES")] <- median(listES)
#   metaDataFrame[i, c("Enrichment.Scores")] <- I(as.list(listES))
#   metaDataFrame[i, c("Gene.List")] <- I(as.list(geneList))
  upperBoundIndex = upperBoundIndex + idIncrements
}

requiredMetaDataFrame <- allDavidIterations[,c(1,3,6,7)]
moltenDataFrame <- melt(data=requiredMetaDataFrame, id=c("Number.Of.Genes"))

ggplot(data=moltenDataFrame, aes(x=moltenDataFrame$Number.Of.Genes, 
                                 y=moltenDataFrame$value, 
                                 group=moltenDataFrame$variable, 
                                 color=moltenDataFrame$variable, 
                                 shape=moltenDataFrame$variable)) + 
  geom_line() + 
  geom_point() + 
  xlab("# Genes") +
  ggtitle("Cluster Characteristics") + guides(color=FALSE)


ggplot(data=moltenDataFrame, aes(x=moltenDataFrame$Number.Of.Genes, 
                                 y=moltenDataFrame$Number.Of.Clusters, 
                                 group=1, 
                                 color="#D55E00")) + 
geom_line() + 
geom_point() +
xlab("# Genes") + ylab("# Cluster") +
ggtitle("Number of Clusters") + guides(color=FALSE)







