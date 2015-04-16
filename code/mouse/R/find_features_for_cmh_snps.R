#!/usr/bin/env Rscript


# find SNPS in Gene Features
# Take in Two aguments directory of
#   - SNP split chromosome wise with columns snp_id, loc
#   - GENE reference csv file
# Assumes that Gene reference file is processed and is a csv file with gene id and corresponding tables
# all paths are absolute
# example call ./

rm(list=ls())
library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

norm_var_names <- function(vars, sep="_") {
  if (sep == ".") sep <- "\\."
  
  # Replace all _ and . and ' ' with the nominated separator.
  
  pat  <- '_|\\.| |,'
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any all capitals words with Initial capitals
  
  pat  <- stringr::perl('(?<!\\p{Lu})(\\p{Lu})(\\p{Lu}*)')
  rep  <- '\\1\\L\\2'
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any capitals not at the beginning of the string with _ 
  # and then the lowercase letter.
  
  pat  <- stringr::perl('(?<!^)(\\p{Lu})')
  rep  <- paste0(sep, '\\L\\1')
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # WHY DO THIS? Replace any number sequences not preceded by an
  # underscore, with it preceded by an underscore. The (?<!...) is a
  # lookbehind operator.
  
  pat  <- stringr::perl(paste0('(?<![', sep, '\\p{N}])(\\p{N}+)'))
  rep  <- paste0(sep, '\\1')
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Remove any resulting initial or trailing underscore or multiples:
  #
  # _2level -> 2level
  
  vars <- stringr::str_replace(vars, "^_+", "")
  vars <- stringr::str_replace(vars, "_+$", "")
  vars <- stringr::str_replace(vars, "__+", "_")
  
  # Convert to lowercase
  
  vars <- tolower(vars)
  
  # Remove repeated separators.
  
  pat  <- paste0(sep, "+")
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  return(vars)
}

# ref_gene_id_data <- NULL
# genome_data_filename <- NULL 
# cmd_arguments <- commandArgs(trailingOnly=TRUE)
# halfWindowSize <- 500000
# ref_gene_id_filename <- NULL
# output_dir_path <- NULL
# 
# 
# for (i in 1:length(cmd_arguments)) {
#     print(paste("arg",as.character(i),"=",cmd_arguments[i]))
#     output_dir_path = cmd_arguments[3]
#     ref_gene_id_filename <- cmd_arguments[2]
#     genome_data_filename <- cmd_arguments[1] 
# }

ref_gene_id_data <- NULL
genome_data_filename <- "/home/data/nehil_multiple_testing_csv/p_adjusted_genome.gwas"
window_size <- 1000
ref_gene_id_filename <- "/home/data/reference/77/Mus_musculus_genes_(GRCm38.p3).csv"
output_dir_path <- "/home/data/"

# 
# system.time(genome_data <- read.csv(genome_data_filename,header=T))
# names(genome_data) <-c("chr_no","snp_pos","ref","p_values","chr_p_adjusted","genome_p_adjusted")
# 
# ref_gene_id_data <- read.csv(file = ref_gene_id_filename, header = TRUE)
# names(ref_gene_id_data) <- norm_var_names(names(ref_gene_id_data))
# ref_gene_id_data$gene_mid_loc = trunc((ref_gene_id_data$gene_start_bp + ref_gene_id_data$gene_end_bp)/2)
# # ref_gene_id_data = ref_gene_id_data[,c("ensembl_gene_id", "chromosome_name", "gene_start_bp", "gene_end_bp", "gene_mid_loc")]
# ref_gene_id_data <- ref_gene_id_data[ , -which(names(ref_gene_id_data) %in% c("description"))]
# chromosome_list <- c(1:19, "X", "Y")
# 
# 
# # create a function to get the feature file name, the gwas file name, a string suffix to be put in outfile names
# 
# 
# 
# find_snps_in_genes_cmh <- function(snp_data, ref_data, output_dir_Path) {
#     
#    
#     cat("Dimesnions of Input Data frames", "Feature Dataset", dim(ref_data), "Snps Dataset", dim(snp_data), "\n\n")    
#     
#     out_file_name <- paste(output_dir_path,"cmh_snps_in_gene_features",".csv", sep="")
#     # Clearing the old version of the filess
#     cat("out_file_name",out_file_name )
#     if (file.exists(out_file_name) == TRUE) {
#         file.remove(out_file_name)
#     }
# 
#     find_store_related_features <- function(snp_row, ref_data, out_file_name) {
#         ref_data$chromosome_name <- as.character(ref_data$chromosome_name)
#         snp_row$chr_no <- as.character(snp_row$chr_no)
#         ref_data$chromosome_name <- str_join("chr", ref_data$chromosome_name)
#         filtered_ref_data <- filter(ref_data, chromosome_name ==  snp_row$chr_no)
#         print(snp_row$chr_co)
#         print(snp_row$snp_pos)
#         cat("Dimesnions::", "Feature Dataset", dim(filtered_ref_data), "Snp Row Chromosome Number", snp_row$chr_no, "\n\n") 
#         features_found_df <- filter(filtered_ref_data,
#                                     (gene_start_bp - window_size) <= snp_row$snp_pos & (gene_end_bp + window_size) >= snp_row$snp_pos
#                                         )
#         print(dim(features_found_df))
#         
#         if (dim(features_found_df)[[1]] > 0) {
#             features_found_df$snp_p_values <- snp_row$p_values
#             features_found_df$snp_pos <- snp_row$snp_pos
#             features_found_df$snp_ref <- snp_row$ref
#             features_found_df$snp_genome_p_values <- snp_row$genome_p_adjusted
#             features_found_df$snp_chr_p_values <- snp_row$chr_p_adjusted
#             print(str(features_found_df))
#             if (file.exists(out_file_name) == TRUE) {
#                 write.table(features_found_df, file = out_file_name, quote=F, sep = ",", append = T,  col.names = FALSE, row.names = FALSE)
#             }  else {
#                 write.table(features_found_df, file = out_file_name, quote=F, sep = ",", col.names = TRUE, row.names = FALSE)
#             }
#             
#         }
#     }
# 
#     ## Big Analysis
#     d_ply(snp_data, .(snp_pos), function(x) {
#         cat("Dimesnions of SNP row", dim(x), "\n\n")
#         find_store_related_features(snp_row = x, 
#                                     ref_data,
#                                     out_file_name)
#     })
# }
# 
# 
# 
# 
# d_ply(genome_data[1:1000,], .(chr_no), function(df) {
#   find_snps_in_genes_cmh(df, ref_gene_id_data, output_dir_path)
# })
# 
# print_dim <- function(dt) {
#   dim(dt)
# }


# The Data Table way
library(data.table)
system.time(genome_data <- fread(genome_data_filename, sep=",", sep2="auto", header=T, na.strings="NA",
      stringsAsFactors=FALSE, verbose=TRUE))
setnames(genome_data,names(genome_data),c("chromosome_name","snp_pos","ref","p_values","chr_p_adjusted","genome_p_adjusted"))
genome_data$fake_gene_start <- genome_data$snp_pos
genome_data$fake_gene_end <- genome_data$snp_pos
ref_gene_id_data <- fread(ref_gene_id_filename, sep=",", sep2="auto", header=T, na.strings="NA",
                        stringsAsFactors=FALSE, verbose=TRUE)
setnames(ref_gene_id_data,names(ref_gene_id_data), norm_var_names(names(ref_gene_id_data)))
ref_gene_id_data$chromosome_name <- str_join("chr", ref_gene_id_data$chromosome_name)
setnames(ref_gene_id_data, c("gene_start_(bp)", "gene_end_(bp)"), c("gene_start", "gene_end"))
ref_gene_id_data$fake_gene_start <- ref_gene_id_data$gene_start - window_size
ref_gene_id_data$fake_gene_end <- ref_gene_id_data$gene_end + window_size

setkey(ref_gene_id_data, chromosome_name, fake_gene_start, fake_gene_end)
result_dt <- foverlaps(genome_data, ref_gene_id_data, type="within", nomatch = 0L)

prod_ready_dt <- result_dt[, .(chromosome_name, ensembl2 * _gene_id, snp_pos, gene_start, gene_end, associated_gene_name, gene_type, ref, p_values, chr_p_adjusted, genome_p_adjusted) ]
prod_ready_dt$ref <- as.character(toupper(as.character(prod_ready_dt$ref)))
write.csv(prod_ready_dt, "snp_in_gene_windows.csv")


df <- fread("snp_in_gene_windows.csv", sep=",", sep2="", header=T)
g_significant <- df[genome_p_adjusted < 0.05]
g_significant <- df[, nlp := -1 * log(genome_p_adjusted)]
qplot(data = g_significant, x = genome_p_adjusted, geom="histogram" )

sync_df <- fread("/home/data/filtered_chr_19.sync", sep="\t", sep2="", header=F)
g_significant <- g_significant[chromosome_name == "chr19"]

setnames(sync_df,names(sync_df),c("chromosome_name","snp_pos","ref","C1M","C2M","S1M","C2F","S2F","S1F","S2M","C1F"))
setkey(sync_df, chromosome_name, snp_pos)

merge2 <- merge(x = g_significant, y = sync_df, by = "snp_pos", all.y = TRUE)
# merge2 <- merge2[1:2000,]
merge2 <- data.table(merge2)
merge2[,(2:8) := NULL]
merge2[,(2) := NULL]
merge2[,(3:4) := NULL]
non_sig_merge2 <- merge2[is.na()]
sig_merge2 <- merge2[!is.na(genome_p_adjusted)]
write.table(non_sig_merge2, "non_significant_chr19.sync", sep="\t", row.names=F, quote = F)
write.table(merge2[375:2000], "sign_and_non_sign_chr19.sync", sep="\t", row.names=F, quote = F)

test <- merge2
split_df <- str_split(merge2$C1M,":")
require(devtools)
source_gist(4676064)

// df - dataframe to be split and merged into (not data table)
// colnames - list of columns to be split
// sep=":"
split_on_character <- function(df, col_name, sep=":") {
  split_df <- str_split(df[,(col_name)],":")
  split_df <- as.data.frame(split_df)
  split_df <- as.data.frame(t(split_df))
  row.names(split_df) <- NULL
  return(split_df)
}

col_names <- c("C1M","C2M","C1F","C2F","S1F","S2M","S1M","S2F")


input_df <- merge2
new_df <- input_df
for (x1 in col_names) {
  out_cols <- split_on_character(as.data.frame(input_df), x1, ":")
  print(dim(out_cols))
  for (x2 in 1:dim(out_cols)[2]) {
    col_name <- paste(x1,x2,sep="_")
    print(col_name)
    new_df[,col_name] <- data.frame(out_cols[,paste("V",x2,sep="")])
  }
}


write.table(new_df, "split_chr19.sync", sep="\t", row.names=F, quote = F)


library(Hmisc)
#implement paper
df <- fread("snp_in_gene_windows.csv", sep=",", sep2="", header=T)
g_significant <- df[genome_p_adjusted < 0.05]
g_significant <- na.omit(g_significant)
g_significant[, nlp := -1 * log(genome_p_adjusted), by="ensembl_gene_id"]
g_significant <- na.omit(g_significant)
g_significant[, mean_nlp := mean(nlp), by="ensembl_gene_id"]
g_significant[, max_nlp := mean(nlp), by="ensembl_gene_id"]

g_significant$quartile <- by(g_significant, g_significant$ensembl_gene_id,  function(x) {cut(x$nlp, 
                                breaks=quantile(x$nlp, probs=seq(0,1, by=0.25)), 
                                include.lowest=TRUE)})
g_significant$quartile <- as.factor(g_significant$quartile)

custom_function <- function(x) {
    return(cut(x, breaks=unique(quantile(x, probs=seq(0,1, by=0.25))), include.lowest=TRUE))
}

cut(c$nlp, 
          breaks=unique(quantile(c$nlp, probs=seq(0,1, by=0.25))), 
          include.lowest=TRUE)
custom_function(c$nlp)

c <- g_significant[ensembl_gene_id == "ENSMUSG00000051951",]

g_significant[, check_length1 := custom_function(nlp), by="ensembl_gene_id"]



ddply(g_significant, )