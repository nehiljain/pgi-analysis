

epigeneticsData <- read.table("Data/RawData/Mouse_Epigenetic_Jan21_2014_AnnotatedFeatures.gff", sep="\t", comment.char="")

ReqColNames <- c('Chromosome.Name', 'Source.Name', 
              'Feature', 'Start.Location','End.Location',
              'Atrributes','Site.Mid.Location')

nameVars <- names(epigeneticsData) %in% c("V6", "V7", "V8")
epigeneticsData <- epigeneticsData[!nameVars]

epigeneticsData["Site.Mid.Location"] <- (epigeneticsData$V4 + epigeneticsData$V5)/2

colnames(epigeneticsData) <- ReqColNames

orderOfCol <- c('Chromosome.Name', 'Source.Name', 
                'Feature', 'Start.Location','End.Location',
                'Site.Mid.Location','Atrributes')

epigeneticsData <- epigeneticsData[, orderOfCol]

epigeneticsData$Chromosome.Name <- as.character(epigeneticsData$Chromosome.Name)
epigeneticsData$Chromosome.Name <- lapply(epigeneticsData$Chromosome.Name, function(x) {strsplit(x, 'r')[[1]][2]})
epigeneticsData$Chromosome.Name <- unlist(epigeneticsData$Chromosome.Name)

write.csv(epigeneticsData, file="Data/ProcessedData/processed_epigenetics_raw_data.csv")


lhevData <- read.csv("Data/ProcessedData/lhEv_result.csv")
lhevData <- merge(lhevData, epigeneticsData, by="Site.Mid.Location", all.x=T)  
write.csv(lhevData,file="Data/ProcessedData/lhEv_result.csv")

testEpiData <- epigeneticsData[sample(100000:1000000,10,rep=TRUE),]
