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

# Create tables
conn <- dbConnect(drv = driver, 'gtex_v7_' %&% population %&% '.db')
dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbExecute(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

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
