#! /usr/bin/env Rscript 

argv <- commandArgs(trailingOnly = TRUE)

infile <- argv[1]
outfile <- argv[2]
outfile_row <- argv[3]

# load modules
suppressWarnings(suppressMessages(library(tidyverse)))

# Read in the data
gene_exp = read.table(file = infile, header = TRUE, sep = "\t" )

# Drop columns we don't need in the gene expression dataframe
gene_exp = gene_exp[-c(1, 3, 4)]

# Rename column one
gene_exp = rename(gene_exp, 'NAME' = Gene_Symbol)

# transpose the gene expression matrix
n = gene_exp$NAME
gene_exp_transpose <- as.data.frame(t(gene_exp[,-1]))
colnames(gene_exp_transpose) <- n

# write out the transposed file in csv format
write.table(gene_exp_transpose, file = outfile, sep = ",", col.names = TRUE, row.names = FALSE)
# write with row names
write.table(gene_exp_transpose, file = outfile_row, sep = ",", col.names = TRUE, row.names = TRUE)
