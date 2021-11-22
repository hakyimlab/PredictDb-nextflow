#! /usr/bin/env Rscript 

argv <- commandArgs(trailingOnly = TRUE)

# Read in the data
peer <- argv[1]
gene_expr <- argv[2]
covs_outfile <- argv[3]
expr_outfile <- argv[4]
resid_outfile <- argv[5]

# Read in the data
peer_factors <- read.csv(file = peer, header = FALSE)
gene_exp_transpose <- read.csv(file = gene_expr, header = TRUE)
#Set the column names for the PEER factors (covariates) as the subject IDs
colnames(peer_factors) <- rownames(gene_exp_transpose)

# write the covariates as a .txt file
write.table(peer_factors, file = covs_outfile, sep = "\t",
            row.names = TRUE)

# linear regression
# copy the matrix
expression <- gene_exp_transpose

#Run linear regression on the PEER factor covariates. 
#Then it sets the residuals to the new expression for that gene.
for (i in 1:length(colnames(gene_exp_transpose))) {
  fit = lm(gene_exp_transpose[,i] ~ t(as.matrix(peer_factors)))
  expression[,i] <- fit$residuals
}

# Write out the final expression file
write.table(expression, file = resid_outfile, sep = "\t",
            row.names = TRUE)

write.table(gene_exp_transpose, file = expr_outfile, sep = "\t",
            row.names = TRUE)
