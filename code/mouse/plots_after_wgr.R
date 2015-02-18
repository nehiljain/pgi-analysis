library(ggplot2)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Function to export plot as pdf, eps and png

ExportPlot <- function(gplot, filename, width=11.69, height=8.27) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
  postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
  print(gplot)
  dev.off()
  png(file = paste(filename, '_.png', sep=""), width = width * 100, height = height * 100)
  print(gplot)
  dev.off()
}

# FileNames

inputClusterGeneFreqFilename = "ProcessedData/top_cluster_stats.csv" 
inputTopGeneListFilename = "ProcessedData/wgr_annotated_top_rank_genelist.csv"


#Input

topClusterGenesFreq <- read.csv("ProcessedData/top_cluster_stats.csv")
topGeneList <- read.csv("FinalData/wgr_annotated_top_rank_genelist.csv")
functionalPathwayResults <- read.table("david_functional_annotation_chart_report.txt",sep="\t")
wholeLhSnpList <-  read.csv('ProcessedData/merged_lhsnp_data_from_alex_files.csv',
                            header=TRUE, comment.char="", nrows=44000)


WholeGeneList <- read.csv( file="ProcessedData/no_dup_gene_id_whole_mouse_genome.csv")

plotBiotypeWholeGenome <-ggplot(data = noDupWholeGenome, 
       aes(x=noDupWholeGenome$Gene.Biotype,
           fill=noDupWholeGenome$Gene.Biotype)) + geom_bar(colour="black", stat="bin") + 
  xlab("Biotype") + ylab("Total No. of Genes") +
  ggtitle("Biotype Gene Density Plot(Whole Genome)")




topClusterGenesFreq$Cluster.Index <- as.factor(topClusterGenesFreq$Cluster.Index)

plotClusterGeneUnique <- ggplot(data = topClusterGenesFreq, 
       aes(x=topClusterGenesFreq$Cluster.Index, 
           y=topClusterGenesFreq$Number.of.UniqueGenes,
           fill=topClusterGenesFreq$Cluster.Index)) + geom_bar(colour="black", stat="identity") + 
  xlab("Cluster Number") + ylab("No. of Unique Genes") +
  ggtitle("Cluster Gene Density Plot")



plotClusterGeneNonUnique <- ggplot(data = topClusterGenesFreq, 
       aes(x=topClusterGenesFreq$Cluster.Index, 
           y=topClusterGenesFreq$Total.Number.of.Genes,
           fill=topClusterGenesFreq$Cluster.Index)) + geom_bar(colour="black", stat="identity") + 
  xlab("Cluster Number") + ylab("Total No. of Genes(non-unique)") +
  ggtitle("Cluster Total Gene Density Plot")


plotBiotypeTotalGenes <- ggplot(data = topGeneList, 
       aes(x=topGeneList$Gene.Biotype,
           fill=topGeneList$Gene.Biotype)) + geom_bar(colour="black", stat="bin") + 
  xlab("Gene Biotype") + ylab("Total No. of Top Genes") +
  ggtitle("Biotype Top-Rank Gene Density Plot")
# 
# functionalAnnotationChartResult <- functionalAnnotationChartResult[with(functionalAnnotationChartResult,order(functionalAnnotationChartResult$Count)),]
# ggplot(data = functionalAnnotationChartResult, 
#        aes(x=functionalAnnotationChartResult$Term,
#            y=functionalAnnotationChartResult$Count)) + geom_point(shape=1) + 
#   xlab("Biological Pathway") + ylab("Total No. of Genes") +
#   ggtitle("Biological Pathway Gene Density Plot")



plotLHGV <- ggplot(wholeGeneList, aes(x=wholeGeneList$Chromosome.Name, y=wholeGeneList$LHGV, 
                                      color=factor(wholeGeneList$Chromosome.Name))) +
  geom_point(shape=19) +
  xlab("Chromosome Number") + ylab("LHGV Values") +
  ggtitle("LHGV Chromosomewise Distributed") + guides(color=FALSE)

plotWGR <- ggplot(wholeGeneList, aes(x=wholeGeneList$Chromosome.Name, y=wholeGeneList$WGR.Score, 
                                     color=factor(wholeGeneList$Chromosome.Name))) +
  geom_point(shape=19) +
  xlab("Chromosome Number") + ylab("WGR Values") +
  ggtitle("WGR Chromosomewise Distributed") + guides(color=FALSE)

plotLHSNP <- ggplot(wholeLhSnpList, aes(x=wholeLhSnpList$CHR, y=wholeLhSnpList$lh, 
                                        color=factor(wholeLhSnpList$CHR))) +
  geom_point(shape=19) +
  xlab("Chromosome Number") + ylab("lhSNP Values") +
  ggtitle("lhSNP Chromosomewise Distributed") + guides(color=FALSE)

multiplot(plotLHSNP, plotLHGV, plotWGR, cols=1)




plotRegWgrLHGV <- ggplot(topGeneList, aes(x=topGeneList$LHGV, y=topGeneList$WGR.Score)) +
  geom_point(shape=1)+
  geom_smooth(method="lm")


