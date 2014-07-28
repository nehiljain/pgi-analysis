Gene Meta Analysis
========================================================

This is an R Markdown document. This is to document all the analysis done and the explaination of the attributes of ach data table created.

Input:

Snps found in features. 
File Name : "SNPS_IN_FEATURES.csv"


```r
library(ggplot2)
snpsInFeatures <- read.csv("~/coderepo/pgi-ngs-analysis/data/final/SNPS_IN_FEATURES.csv", 
    nrows = 10, stringsAsFactors = F, row.names = NULL, header = T)

head(snpsInFeatures)
```

```
##   row.names chromosome.no               source    feature   start     end
## 1         1          chr1       protein_coding       gene 3205901 3671498
## 2         2          chr1 processed_transcript transcript 3205901 3216344
## 3         3          chr1 processed_transcript       exon 3205901 3207317
## 4         1          chr1       protein_coding       gene 3205901 3671498
## 5         2          chr1 processed_transcript transcript 3205901 3216344
## 6         3          chr1 processed_transcript       exon 3205901 3207317
##   score strand frame
## 1     .      -     .
## 2     .      -     .
## 3     .      -     .
## 4     .      -     .
## 5     .      -     .
## 6     .      -     .
##                                                                                                                                                                                                                                                                                                    attributes
## 1                                                                                                                                                                                                        gene_id ENSMUSG00000051951; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding;
## 2                                            gene_id ENSMUSG00000051951; transcript_id ENSMUST00000162897; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Xkr4-003; transcript_source havana; tag cds_end_NF; tag cds_start_NF; tag mRNA_end_NF; tag mRNA_start_NF;
## 3 gene_id ENSMUSG00000051951; transcript_id ENSMUST00000162897; exon_number 2; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Xkr4-003; transcript_source havana; exon_id ENSMUSE00000866652; tag cds_end_NF; tag cds_start_NF; tag mRNA_end_NF; tag mRNA_start_NF;
## 4                                                                                                                                                                                                        gene_id ENSMUSG00000051951; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding;
## 5                                            gene_id ENSMUSG00000051951; transcript_id ENSMUST00000162897; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Xkr4-003; transcript_source havana; tag cds_end_NF; tag cds_start_NF; tag mRNA_end_NF; tag mRNA_start_NF;
## 6 gene_id ENSMUSG00000051951; transcript_id ENSMUST00000162897; exon_number 2; gene_name Xkr4; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Xkr4-003; transcript_source havana; exon_id ENSMUSE00000866652; tag cds_end_NF; tag cds_start_NF; tag mRNA_end_NF; tag mRNA_start_NF;
##                      gene_id      gene_name                gene_biotype
## 1 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
## 2 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
## 3 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
## 4 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
## 5 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
## 6 gene_id ENSMUSG00000051951 gene_name Xkr4 gene_biotype protein_coding
##   snp_name  snp_bp snp_p_value
## 1   snp956 3206270   1.827e-04
## 2   snp956 3206270   1.827e-04
## 3   snp956 3206270   1.827e-04
## 4   snp957 3206491   1.040e-46
## 5   snp957 3206491   1.040e-46
## 6   snp957 3206491   1.040e-46
```



