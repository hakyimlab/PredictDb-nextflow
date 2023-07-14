#! /usr/bin/env Rscript 

## Process the peer covariates
argv <- commandArgs(trailingOnly = TRUE)

# Read in the data
comp_cov <- argv[1]
prov_cov <- argv[2]
covs_outfile <- argv[3]

# Read in the data
comp_cov <- read.table(comp_cov, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
prov_cov <- read.table(prov_cov, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#find the common columns
common_cols <- intersect(colnames(comp_cov), colnames(prov_cov))

# order the columns
comp_cov <- comp_cov[,common_cols]
prov_cov <- prov_cov[,common_cols]

# bind the covariate data
all_covs <- rbind(comp_cov, prov_cov)

# write the covariates as a .txt file
write.table(all_covs, file = covs_outfile, sep = "\t",
            row.names = TRUE)

