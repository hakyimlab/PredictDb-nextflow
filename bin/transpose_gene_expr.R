#! /usr/bin/env Rscript 

argv <- commandArgs(trailingOnly = TRUE)

infile <- argv[1]
outfile <- argv[2]
outfile_tsv <- argv[3]

# load modules
suppressWarnings(suppressMessages(library(tidyverse)))

# Read in the data in a flexible manner
if (grepl("\\.txt$", infile)) {
  gene_exp = read.table(file = infile, header = TRUE, sep = "\t" )
} else if (grepl("\\.tsv$", infile)) {
  gene_exp = read.table(file = infile, header = TRUE, sep = "\t" )
} else if (grepl("\\.csv$", infile)) {
  gene_exp <- read.csv(infile, header = TRUE)
} else {
  stop("Invalid file format, check the extension")
}

# Drop columns we don't need in the gene expression dataframe
drop_cols <- c("Chr","TargetID","Coord")

gene_exp <- gene_exp[,!(names(gene_exp) %in% drop_cols)]

# Rename column one
gene_exp <- gene_exp %>% rename(Gene_name = names(.)[1])

# transpose the gene expression matrix
n = gene_exp$Gene_name
gene_exp_transpose <- as.data.frame(t(gene_exp[,-1]))
colnames(gene_exp_transpose) <- n

# write out the transposed file in csv format for peer factors
write.table(gene_exp_transpose, file = outfile, sep = ",", 
            col.names = TRUE, row.names = TRUE)

# write out the transposed file in tsv format for downstream analysis
write.table(gene_exp_transpose, file = outfile_tsv, sep = "\t",
            row.names = TRUE)
