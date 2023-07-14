#! /usr/bin/env Rscript 

## Process the peer covariates

argv <- commandArgs(trailingOnly = TRUE)

# Read in the data
peer <- argv[1]
gene_expr <- argv[2]
covs_outfile <- argv[3]

# Read in the data
peer_factors <- read.csv(file = peer, header = FALSE)
gene_exp_transpose <- read.csv(file = gene_expr, header = TRUE)
#Set the column names for the PEER factors (covariates) as the subject IDs
colnames(peer_factors) <- rownames(gene_exp_transpose)

# write the covariates as a .txt file
write.table(peer_factors, file = covs_outfile, sep = "\t",
            row.names = TRUE)

