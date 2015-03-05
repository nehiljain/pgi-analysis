#Merge LhSnp Values from Alex Files


lh_X <- read.table("../..//Data//RawData//alex_lhgv_calculation_files/linkage/lh_X.txt", sep="\t", comment.char="")
lh_Y <- read.table("../..//Data//RawData//alex_lhgv_calculation_files/linkage/lh_Y.txt", sep="\t", comment.char="")
lh.names <- paste("../..//Data//RawData//alex_lhgv_calculation_files/linkage/lh_", 1:19,".txt", sep = "")

lh_1 <- read.table(lh.names[1], sep="\t", comment.char="")
lh_2 <- read.table(lh.names[2], sep="\t", comment.char="")
lh_3 <- read.table(lh.names[3], sep="\t", comment.char="")
lh_4 <- read.table(lh.names[4], sep="\t", comment.char="")
lh_5 <- read.table(lh.names[5], sep="\t", comment.char="")
lh_6 <- read.table(lh.names[6], sep="\t", comment.char="")
lh_7 <- read.table(lh.names[7], sep="\t", comment.char="")
lh_8 <- read.table(lh.names[8], sep="\t", comment.char="")
lh_9<- read.table(lh.names[9], sep="\t", comment.char="")
lh_10 <- read.table(lh.names[10], sep="\t", comment.char="")
lh_11 <- read.table(lh.names[11], sep="\t", comment.char="")
lh_12 <- read.table(lh.names[12], sep="\t", comment.char="")
lh_13 <- read.table(lh.names[13], sep="\t", comment.char="")
lh_14 <- read.table(lh.names[14], sep="\t", comment.char="")
lh_15 <- read.table(lh.names[15], sep="\t", comment.char="")
lh_16 <- read.table(lh.names[16], sep="\t", comment.char="")
lh_17 <- read.table(lh.names[17], sep="\t", comment.char="")
lh_18 <- read.table(lh.names[18], sep="\t", comment.char="")
lh_19 <- read.table(lh.names[19], sep="\t", comment.char="")

lh_1 ["CHR"]<- 1
lh_2 ["CHR"]<- 2
lh_3 ["CHR"]<- 3
lh_4 ["CHR"]<- 4
lh_5 ["CHR"]<- 5
lh_6 ["CHR"]<- 6
lh_7 ["CHR"]<- 7
lh_8 ["CHR"]<- 8
lh_9 ["CHR"]<- 9
lh_10 ["CHR"]<- 10
lh_11 ["CHR"]<- 11
lh_12 ["CHR"]<- 12
lh_13 ["CHR"]<- 13
lh_14 ["CHR"]<- 14
lh_15 ["CHR"]<- 15
lh_16 ["CHR"]<- 16
lh_17 ["CHR"]<- 17
lh_18 ["CHR"]<- 18
lh_19 ["CHR"]<- 19
lh_X ["CHR"] <- 20
lh_Y ["CHR"] <- 21

wholeLhSnpData <- rbind(lh_1 ,
lh_2 ,
lh_3 ,
lh_4 ,
lh_5 ,
lh_6 ,
lh_7 ,
lh_8 ,
lh_9 ,
lh_10 ,
lh_11 ,
lh_12 ,
lh_13 ,
lh_14 ,
lh_15 ,
lh_16 ,
lh_17 ,
lh_18 ,
lh_19 ,

lh_X, lh_Y)


write.csv(wholeLhSnpData, file="ProcessedData/merged_lhsnp_data_from_alex_files.csv",sep=",")

