#  "Wed Jan  1 20:39:55 2014"


wholeGenomeFname = "../..//Data//RawData//Mouse_Genome_Features_Dec_16_2013.csv"
noDuplicateWHoleGenomeOutputFname = "ProcessedData/no_dup_gene_id_whole_mouse_genome.csv"
wholeGenomeforlhSNPCalcOutputFname= "ProcessedData/whole_mouse_genome_for_lhSnp_calc.csv"


initialRows = read.csv(wholeGenomeFname, nrow = 1000)


wholeGenomeClasses <- sapply(initialRows, class)


wholeGenome <- read.csv(wholeGenomeFname, header=TRUE, comment.char="", nrows=900000)

head(wholeGenome)


#REMOVE ALL THE DUPLICATE ROWS DISCARDING THE MANY TRANSCRIPT IDS ETC
# relevantWholeGenome <- relevantWholeGenome[!duplicated(relevantWholeGenome$Ensembl.Gene.ID),]

colRequired <- as.logical(c("TRUE", "FALSE" , "TRUE", "TRUE", "TRUE" ))

# the part of genome relevant for processing of LHGV file
relevantWholeGenome <- wholeGenome[1:6] 


#REMOVE ALL THE DUPLICATE ROWS DISCARDING THE MANY TRANSCRIPT IDS ETC
noDupWholeGenome <- wholeGenome[!duplicated(wholeGenome$Ensembl.Gene.ID),]


noDupWholeGenome$Gene.Biotype <- as.character(noDupWholeGenome$Gene.Biotype)




noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "protein_coding" ] <- "Amino_Acid_Coding" 
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "pseudogene" ] <- "Amino_Acid_Coding"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "polymorphic_pseudogene" ] <- "Amino_Acid_Coding"

noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "antisense" ] <- "Other" 
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "sense_intronic" ] <- "Other"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "sense_overlapping" ] <- "Other"


noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_LV_gene" ] <- "Intergenic"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_D_gene"] <- "Intergenic"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_J_gene"] <- "Intergenic"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_V_gene"] <- "Intergenic"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_C_gene"] <- "Intergenic"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "IG_V_pseudogene" ] <- "Intergenic"

noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "processed_transcript" ] <- "Transcript"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "TR_V_gene" ] <- "Transcript"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "TR_V_pseudogene" ] <- "Transcript"

noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "Mt_tRNA"] <- "Mitochondrial"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "Mt_rRNA" ] <- "Mitochondrial"

noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "lincRNA" ] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "3prime_overlapping_ncrna" ] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "snoRNA" ] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "miRNA" ] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "miscRNA" ] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "snRNA"] <- "RNA"
noDupWholeGenome$Gene.Biotype[noDupWholeGenome$Gene.Biotype == "rRNA" ] <- "RNA"

write.csv(noDupWholeGenome, file=noDuplicateWHoleGenomeOutputFname)


lhSnpWholeGenome <- relevantWholeGenome[colRequired]

#calculate mid location
lhSnpWholeGenome <- within(lhSnpWholeGenome, Mid.Location <- (Gene.End..bp. + Gene.Start..bp.) * 0.5)

lhSnpWholeGenome <- unique(lhSnpWholeGenome)
lhSnpWholeGenome <- lhSnpWholeGenome[with(lhSnpWholeGenome, order(Chromosome.Name)),]
write.csv(lhSnpWholeGenome, file=wholeGenomeforlhSNPCalcOutputFname)

rm(relevantWholeGenome, wholeGenome)

