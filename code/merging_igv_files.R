library(stringr)
setwd("~/Desktop/cmh-igv-files/")
filelist <- list.files("~/Desktop/cmh-igv-files/")

for (filename in filelist){ 
  print(filename)
  varname <- str_split(filename,":", n=2)[[1]][1]
  assign(varname, read.table(filename, header=T, sep="\t", stringsAsFactors=F))
}

merged <- rbind(allSamples.chr1,allSamples.chr10, allSamples.chr11, allSamples.chr12, 
                allSamples.chr13, allSamples.chr14, allSamples.chr15, allSamples.chr16,
                allSamples.chr17, allSamples.chr18, allSamples.chr19, allSamples.chr2,
                allSamples.chr3, allSamples.chr4, allSamples.chr5, allSamples.chr6, 
                allSamples.chr7, allSamples.chr8, allSamples.chr9, allSamples.chrX, allSamples.chrY)
merged[merged$P == 0, 4] <- 0.0000000000000000000001
write.table(merged, "../merged-cmh-igv-all-chromosomes.gwas", sep="\t", row.names=F, quote=F)
