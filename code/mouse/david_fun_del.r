library(RDAVIDWebService)


# Constant PathVariables
clusterReportFname = "david_cluster_report.txt"
funcAnnotTableFname = "david_functional_annotation_table_report.txt"
funcAnnotChartFname = "david_functional_annotation_chart_report.txt"
mergedLHGVOutputFname = "ProcessedData/merged_lhgv_results_whole_mouse_genome.csv" 

# Read Genelist 
wholeGeneList <-  read.csv(mergedLHGVOutputFname, 
                            header=TRUE, comment.char="", nrows=40000)

#reducing size to what is relevant  
wholeGeneList <- wholeGeneList[1:3]
missingValue <- is.na(wholeGeneList$LHGV)
wholeGeneList[missingValue,"LHGV"] <- - 99.99
topLhGVBool <- wholeGeneList$LHGV > 1.5
topGeneList <- wholeGeneList[topLhGVBool,"Ensembl.Gene.ID"]


# Connect to DAVID

david <- DAVIDWebService$new(email = 'nh736060@dal.ca')
connect(david)
idTypes <- getIdTypes(david)

addListResult <- addList(david, topGeneList, idType="ENSEMBL_GENE_ID", listName = "mouseGeneList", listType=c("Gene"))

# Generate all the Report Files
annotationSummaryResult <- getAnnotationSummary(david)

geneListReport <- getGeneListReport(david)

# getFunctionalAnnotationChartFile(david, funcAnnotChartFname)
getClusterReportFile(david, clusterReportFname)

# getFunctionalAnnotationTableFile(david, funcAnnotTableFname)


# functionalAnnotationChartResult <- getFunctionalAnnotationChart(david)
clusterReport <- getClusterReport(david)
  
