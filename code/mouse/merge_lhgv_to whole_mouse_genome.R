#merges results from lhgc_results and reorders the columns

# The line below can be omitted if already done.
# source('~/Dropbox/PGI (1)/WGRanalysis/Code/process_whole_mouse_reference_genome_data.r')

#IMP Might need to change the file name

kAnalysisName = 'lhgv_133'


lhgvFilename = 'ProcessedData/lhgv_result.csv'
noDupWholeGenome <-  read.csv('ProcessedData/no_dup_gene_id_whole_mouse_genome.csv')


lhgvResult <- read.csv(lhgvFilename, header=TRUE, comment.char="", nrows=40000)

#LEFT Outer Join is done
mergeResult <- merge(noDupWholeGenome,lhgvResult, by='Ensembl.Gene.ID', all.x=TRUE)

#reordering

orderOfCol <- c("Ensembl.Gene.ID",  "LHGV", "Chromosome.Name", "Gene.Mid.Location", "Gene.Biotype", "Associated.Gene.Name", 
  "Description", "Phenotype.description", 
  "No.Of.Snps", "Top.Lh.Snp", "Top.Snp.Name" , "Top.Snp.Location",  
  "Ensembl.Transcript.ID", "Gene.Start..bp."   ,  "Gene.End..bp.",    "Strand" ,"Band",
  "Transcript.Start..bp." ,    
  "Transcript.End..bp." ,  "Associated.Gene.DB"  ,"Transcript.count" ,"X..GC.content",                                           
  "Ensembl.Exon.ID"  ,          "Ensembl.Protein.ID"   ,      "Associated.Transcript.Name",
  "Associated.Transcript.DB","Transcript.Biotype" ,"Source" ,                   
  "Status..gene."     , "Status..transcript."   ,          
  "Source.name"  , "Study.External.Reference" ,"Strain.name",               
  "Strain.gender" ,"P.value" )        
    

reOrderedMergeResult <- mergeResult[orderOfCol]


write.csv(reOrderedMergeResult, file="ProcessedData/merged_lhgv_results_whole_mouse_genome.csv")
