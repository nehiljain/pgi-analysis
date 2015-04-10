
df <- fread("snp_in_gene_windows.csv", sep=",", sep2="", header=T)
g_significant <- df[genome_p_adjusted < 0.05]

g_significant[, min_p := min(genome_p_adjusted), by="ensembl_gene_id"]
g_significant[, nlp := -1 * log(genome_p_adjusted), by="ensembl_gene_id"]
g_significant[, max_nlp := max(nlp), by="ensembl_gene_id"]
g_significant[, mean_nlp := mean(nlp), by="ensembl_gene_id"]





calculate_max_nlp_per_gene <- function(df) {
  
}