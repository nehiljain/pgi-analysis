
library(tm)
library(wordcloud)
library(RColorBrewer)
#read the file

topClusterData <- read.csv("ProcessedData/top_cluster_stats.csv")

topClusterData$Biological.Terms <- as.character(topClusterData$Biological.Terms)
topClusterData$Biological.Terms <- strsplit(topClusterData$Biological.Terms, ",")

topClusterData$Biological.Categories <- as.character(topClusterData$Biological.Categories)
topClusterData$Biological.Categories <- strsplit(topClusterData$Biological.Categories , ",")

createWordCloud <- function(index){
  ap.corpus <- Corpus(DataframeSource(data.frame(as.character(topClusterData$Biological.Terms[index]))))
  ap.corpus <- tm_map(ap.corpus, removePunctuation)
  ap.corpus <- tm_map(ap.corpus, tolower)
  
  ap.tdm <- TermDocumentMatrix(ap.corpus)
  ap.m <- as.matrix(ap.tdm)
  ap.v <- sort(rowSums(ap.m),decreasing=TRUE)
  ap.d <- data.frame(word = names(ap.v),freq=ap.v)
  table(ap.d$freq)
  pal2 <- brewer.pal(8,"Dark2")
  wordcloudFname = paste("../Figures/wordcloud_packages_",index,".png", sep="")
  png(wordcloudFname, width=1280,height=800)
  wordcloud(ap.d$word,ap.d$freq, scale=c(4,.5),min.freq=1,
            max.words=Inf, random.order=FALSE, rot.per=.15, colors=pal2)
  dev.off()
}


wordCloudResults <- lapply(topClusterData$Cluster.Index, function(x) createWordCloud(x))
