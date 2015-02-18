
davidGeneList <-  read.csv('ProcessedData/wgr_result.csv')
wholeGeneList <-  read.csv('ProcessedData/merged_lhgv_results_whole_mouse_genome.csv',
                           header=TRUE, comment.char="", nrows=40000)


davidGeneList$CS <- davidGeneList$TempNumeratorCSCalc/davidGeneList$CSFREQ

minLHGV <- min(wholeGeneList$LHGV, na.rm = T)
maxLHGV <- max(wholeGeneList$LHGV, na.rm = T)
rangeLHGV = maxLHGV -minLHGV
wholeGeneList['Normalised.LHGV'] <- (wholeGeneList$LHGV-minLHGV) / rangeLHGV

davidGeneList <- davidGeneList[, c("Ensembl.Gene.ID", 
                                   "Name.of.Included.FunctionalPathways","Name.of.Excluded.FunctionalPathways", "CS","CSFREQ"
                                   ,"FuncFreq", "FP", "Total.Number.Of.Pathways" ) ]

colnames(davidGeneList) <- c("Ensembl.Gene.ID", "Name.of.Included.FunctionalPathways", 
                             "Name.of.Excluded.FunctionalPathways", 
                             "Cluster.Score", "No.Of.Custers", 
                             "No.Of.FunctionalPathways", "FunctionalPathway.Score", "Total.Number.Of.Pathways")
davidGeneList <- davidGeneList[, c("Ensembl.Gene.ID", 
                             "Cluster.Score", "No.Of.Custers", 
                             "FunctionalPathway.Score", "No.Of.FunctionalPathways", "Total.Number.Of.Pathways","Name.of.Included.FunctionalPathways", 
                             "Name.of.Excluded.FunctionalPathways")]
wgrGeneList <- merge(x = wholeGeneList, y = davidGeneList, by = "Ensembl.Gene.ID",all.y = TRUE)

wgrGeneList$Cluster.Score[is.na(wgrGeneList$Cluster.Score)] <- 0
wgrGeneList$FunctionalPathway.Score[is.na(wgrGeneList$FunctionalPathway.Score)] <- 0



stdLHGV = sd(wgrGeneList$LHGV)
stdCS = sd(wgrGeneList$Cluster.Score)
stdFPS = sd(wgrGeneList$FunctionalPathway.Score)
wgrGeneList["WGR.Score"] <- ((0.75 * wgrGeneList$LHGV/stdLHGV) 
                             + (0.15 * wgrGeneList$Cluster.Score/stdCS)
                             + (0.10 * (wgrGeneList$FunctionalPathway.Score/(wgrGeneList$No.Of.FunctionalPathways * stdFPS)))) 

wgrGeneList$WGR.Score[is.na(wgrGeneList$WGR.Score)] <- ((0.75 * wgrGeneList$LHGV/stdLHGV) 
                                                        + (0.15 * wgrGeneList$Cluster.Score/stdCS))
                                                      

write.csv(wgrGeneList, file="FinalData/wgr_annotated_top_rank_genelist.csv")


