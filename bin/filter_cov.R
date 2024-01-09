#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))

args <- commandArgs(trailingOnly = TRUE)

filtered_db <- args[1]
unfiltered_cov <- args[2]
filtered_cov <- args[3]

# load db
driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from extra')
dbDisconnect(in_conn)

# load covariances
all_covs <- read.table(unfiltered_cov, header = TRUE, stringsAsFactors = FALSE)

# Filter out models with low performance
all_covs <- all_covs %>% 
    filter(GENE %in% model_summaries$gene)

# Write out the covariance file filtered
write.table(all_covs, file = filtered_cov, sep = "\t",
            row.names = FALSE, quote = FALSE)
