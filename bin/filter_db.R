#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))

args <- commandArgs(trailingOnly = TRUE)


unfiltered_db <- args[1]
filtered_db <- args[2]

driver <- dbDriver("SQLite")

in_conn <- dbConnect(driver, unfiltered_db)
out_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from extra')

# Filter out models with low performance
model_summaries <- model_summaries %>% 
    filter(pred.perf.pval < 0.05 & rho_avg > 0.1) %>%
    filter(n.snps.in.model > 0) 


dbWriteTable(out_conn, 'extra', model_summaries, overwrite = TRUE)
construction <- dbGetQuery(in_conn, 'select * from construction')
dbWriteTable(out_conn, 'construction', construction, overwrite = TRUE)
sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
dbWriteTable(out_conn, 'sample_info', sample_info, overwrite = TRUE)

weights <- dbGetQuery(in_conn, 'select * from weights')
weights <- weights %>% filter(gene %in% model_summaries$gene)
dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
dbExecute(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbExecute(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
dbExecute(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbExecute(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")

dbDisconnect(in_conn)
dbDisconnect(out_conn)
