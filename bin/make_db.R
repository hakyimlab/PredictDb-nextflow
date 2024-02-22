#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))

"%&%" <- function(a,b) paste(a,b,sep='')

args <- commandArgs(trailingOnly = TRUE)

model_summary <- args[1]
weight_summary <- args[2]
chroms_summary <- args[3]
population <- args[4]

driver <- dbDriver('SQLite')

# Read in the files
model_summaries <- read.table(model_summary, header = T, stringsAsFactors = F)
weights <- read.table(weight_summary, header = T, stringsAsFactors = F)
chrom_summary <- read.table(chroms_summary, header = T, stringsAsFactors = F)

# Rename columns
if('rho_avg_squared' %in% colnames(model_summaries)){
    # nested cv elasticnet
    model_summaries <- model_summaries %>%
        rename(pred.perf.R2 = rho_avg_squared, gene = gene_id,
            genename = gene_name, pred.perf.pval = zscore_pval, 
            n.snps.in.model = n_snps_in_model)%>%
        dplyr::mutate(pred.perf.qval = NA)

} else {
    # elasticnet
    model_summaries <- model_summaries %>%
        rename(pred.perf.R2 = rho_avg_squared, gene = gene_id,
            genename = gene_name, pred.perf.pval = zscore_pval, 
            n.snps.in.model = n_snps_in_model) %>%
        dplyr::mutate(pred.perf.qval = NA)
}

weights <- weights %>%
    rename(eff_allele = alt, ref_allele = ref, weight = beta, gene = gene_id)

# Create tables
conn <- dbConnect(drv = driver, 'predict_db_' %&% population %&% '.db')
dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON extra (gene)")

dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
n_samples <- chrom_summary$n_samples
sample_info <- data.frame(n_samples = n_samples, population = population)
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
construction <- chrom_summary %>%
                  select(chrom, cv_seed) %>%
                  rename(chromosome = chrom)
dbWriteTable(conn, 'construction', construction, overwrite = TRUE)

dbDisconnect(conn)
