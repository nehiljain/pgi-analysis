require("ggplot2")
require(profr)
require(ff)

chr10 <- read.table.ffdf(x = NULL, "data/popoolation/allSamples.chr10:1-130694993.fst", VERBOSE=T, first.rows=100)

#Getting list of FST Files

#List of all files
allFileList <- list.files(path="data/", recursive = T)

#Logical vector to filter FST files
fstFileList <- grepl(pattern = "fst$", allFileList,ignore.case = T )

#List of FST Files to read in 
fstFileList <- allFileList[fstFileList]

#for each file 

for(file in fstFileList) {
  #getting chromosome number  
  print(file)
  startPos <- regexpr(pattern = "chr", file, ignore.case=T)
  endPos <- regexpr(pattern = ":", file, ignore.case=T)
  
  print(substr(file, startPos[1],endPos[1]-1))

}

