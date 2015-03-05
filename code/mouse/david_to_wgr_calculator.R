rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

davidGeneList <-  read.csv('/home/data/microarray_80k/nehil_pathway_david/wgr_result.csv', stringsAsFactors = FALSE)
wholeGeneList <-  read.csv('/home/data/microarray_80k/nehil_pathway_david/pathway-genes.csv', header=TRUE)

davidGeneList$CS = as.numeric(davidGeneList$CS)
davidGeneList$TempNumeratorCSCalc = as.numeric(davidGeneList$TempNumeratorCSCalc)
davidGeneList$CSFREQ = as.numeric(davidGeneList$CSFREQ)
davidGeneList$CS <- davidGeneList$TempNumeratorCSCalc/davidGeneList$CSFREQ


davidGeneList <- davidGeneList[, c("Ensembl.Gene.ID", 
                                   "Name.of.Included.FunctionalPathways","Name.of.Excluded.FunctionalPathways", "CS","CSFREQ"
                                   ,"FuncFreq", "FP", "Total.Number.Of.Pathways" ) ]

colnames(davidGeneList) <- c("Ensembl.Gene.ID", "Name.of.Included.FunctionalPathways", 
                             "Name.of.Excluded.FunctionalPathways", 
                             "Cluster.Score", "No.Of.Custers", 
                             "No.Of.FunctionalPathways", "FunctionalPathway.Score", "Total.Number.Of.Pathways")
names(davidGeneList) =  c("ensembl_gene_id", "Name.of.Included.FunctionalPathways", 
                          "Name.of.Excluded.FunctionalPathways", 
                          "Cluster.Score", "No.Of.Custers", 
                          "No.Of.FunctionalPathways", "FunctionalPathway.Score", "Total.Number.Of.Pathways")
davidGeneList <- davidGeneList[, c("ensembl_gene_id", 
                             "Cluster.Score", "No.Of.Custers", 
                             "FunctionalPathway.Score", "No.Of.FunctionalPathways", "Total.Number.Of.Pathways","Name.of.Included.FunctionalPathways", 
                             "Name.of.Excluded.FunctionalPathways")]

wgrGeneList <- merge(x = wholeGeneList, y = davidGeneList, by = "ensembl_gene_id",all.y = TRUE)

wgrGeneList$Cluster.Score[is.na(wgrGeneList$Cluster.Score)] <- 0
wgrGeneList$FunctionalPathway.Score[is.na(wgrGeneList$FunctionalPathway.Score)] <- 0


minlhgv <- min(wgrGeneList$lhgv, na.rm = T)
maxlhgv <- max(wgrGeneList$lhgv, na.rm = T)
rangelhgv = maxlhgv -minlhgv
wgrGeneList['norm_lhgv'] =  (wgrGeneList$lhgv-minlhgv) / rangelhgv

min_cs <- min(wgrGeneList$Cluster.Score)
max_cs <- max(wgrGeneList$Cluster.Score)
range_cs = max_cs -min_cs
wgrGeneList['norm_cs'] =  (wgrGeneList$Cluster.Score-min_cs) / range_cs

min_fps <- min(wgrGeneList$FunctionalPathway.Score)
max_fps <- max(wgrGeneList$FunctionalPathway.Score)
range_fps = max_fps -min_fps
wgrGeneList['norm_fps'] =  (wgrGeneList$FunctionalPathway.Score-min_fps) / range_fps

# 
# 
mean_lhgv = mean(wgrGeneList$lhgv)
mean_cs = mean(wgrGeneList$Cluster.Score)
mean_fps = mean(wgrGeneList$FunctionalPathway.Score)

sdlhgv = sd(wgrGeneList$lhgv)
sdcs = sd(wgrGeneList$Cluster.Score)
sdfps = sd(wgrGeneList$FunctionalPathway.Score)
 
wgrGeneList['std_lhgv'] =  (wgrGeneList$lhgv-mean_lhgv) / sdlhgv
wgrGeneList['std_cs'] =  (wgrGeneList$Cluster.Score-mean_cs) / sdcs
wgrGeneList['std_fps'] =  (wgrGeneList$FunctionalPathway.Score) / sdfps


result <- ddply(wgrGeneList, .(ensembl_gene_id), function(row) {
  if (row$Cluster.Score == 0) {
    row$WGR.Score <- (0.75 * (row$lhgv - mean_lhgv)/sdlhgv)
#     print(row)
  } else if (row$No.Of.FunctionalPathways == 0) {
    row$WGR.Score <- (0.75 * (row$lhgv - mean_lhgv)/sdlhgv) + (0.15 * (row$Cluster.Score - mean_cs)/sdcs)
#     print(row)                    
  } else if (row$No.Of.FunctionalPathways > 0 && row$Cluster.Score > 0 ) {
    row$WGR.Score <- (
                      (0.75 * (row$lhgv - mean_lhgv)/sdlhgv) + 
                      (0.15 * (row$Cluster.Score - mean_cs)/sdcs) + 
                      (0.10 * (row$FunctionalPathway.Score-mean_fps)/sdfps)
                      )
#     print(row)
  } else {
#     print(row)
  }
  return(row)
})

min_wgr = min(result$WGR.Score)
max_wgr = max(result$WGR.Score)
range_wgr = max_wgr - min_wgr
result['norm_wgr'] =  (result$WGR.Score-min_wgr) / range_wgr
result = arrange(result, desc(WGR.Score))

result[1:100, c('rank_label')] = rep('top',100)   
result[101:1977, c('rank_label')] = rep('bottom', 1877)

p1 = ggplot(result, aes(x=WGR.Score, y=std_lhgv, color=factor(rank_label))) + geom_point() + theme(legend.position=c(0,0),legend.justification=c(1,1))

p2 = ggplot(result, aes(x=WGR.Score, y=std_cs, color=factor(rank_label))) + geom_point() + theme(legend.position=c(1,1),legend.justification=c(1,1))

p3 = ggplot(result, aes(x=WGR.Score, y=std_fps, color=factor(rank_label))) + geom_point() + theme(legend.position=c(1,1),legend.justification=c(1,1))

p4 = ggplot(result, aes(x=std_lhgv, y=std_fps, color=factor(rank_label))) + geom_point() + theme(legend.position=c(1,1),legend.justification=c(1,1))

p5 = ggplot(result, aes(x=std_lhgv, y=std_cs, color=factor(rank_label))) + geom_point() + theme(legend.position=c(1,1),legend.justification=c(1,1))

p6 = ggplot(result, aes(x=std_fps, y=std_cs, color=factor(rank_label))) + geom_point() + theme(legend.position=c(1,1),legend.justification=c(1,1))

a <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow=2)
write.csv(result , file="/home/data/microarray_80k/wgr_annotated_top_rank_genelist.csv")


