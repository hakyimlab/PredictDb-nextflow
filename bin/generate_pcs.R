#! /usr/bin/env Rscript 

## Process the peer covariates
argv <- commandArgs(trailingOnly = TRUE)

# Read in the data
gene_expr <- argv[1]
covs_outfile <- argv[2]

# Read in the data
gene_exp_transpose <- read.csv(file = gene_expr, header = TRUE)

# compute the PCs
gene_pcs <- prcomp(gene_exp_transpose, scale. = F, center = T)

# extract the top 10 PCs
pcs <- gene_pcs$x
pcs <- pcs[,1:10]
rownames(pcs) <- rownames(gene_exp_transpose)

# transpose the PCs
pcs <- t(pcs)

# write the covariates as a .txt file
write.table(pcs, file = covs_outfile, sep = "\t",
            row.names = TRUE)

