# File to make manhattan plots for PGI
library(ggplot2)
# par(mfcol = c(3,1))
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
wholeGeneList <-  read.csv("./Data/ProcessedData/wgr_annotated_top_rank_genelist.csv",
                           header=TRUE, comment.char="", nrows=40000)

wholeLhSnpList <-  read.csv('Data/ProcessedData/merged_lhsnp_data_from_alex_files.csv',
                           header=TRUE, comment.char="", nrows=44000)

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
# source('./Code/qqman.r')
# 
# lhgvManhattanData <- wholeGeneList[2:5]
# missingValue <- is.na(lhgvManhattanData$LHGV)
# lhgvManhattanData[missingValue,"LHGV"] <- - 99.99
# colnames(lhgvManhattanData) <- c( "SNP", "BP", "P",  "CHR")
# lhgvManhattanData$CHR <- as.character(lhgvManhattanData$CHR)
# lhgvManhattanData$CHR[lhgvManhattanData$CHR == "X"] <- "20"
# lhgvManhattanData$CHR[lhgvManhattanData$CHR == "Y"] <- "21"
# lhgvManhattanData$CHR <- as.numeric(lhgvManhattanData$CHR)
# complete.cases(lhgvManhattanData$CHR)
# 
# ggplot(lhgvManhattanData, aes(x=lhgvManhattanData$CHR, y=lhgvManhattanData$P, 
#                               color=factor(lhgvManhattanData$CHR))) +
#   geom_point(shape=19) 
# 
# manhattan(lhgvManhattanData, colors=c("black","#666666","#CC6600"), pch=20, genomewideline=F, suggestiveline=F)
# 
# 
# 

