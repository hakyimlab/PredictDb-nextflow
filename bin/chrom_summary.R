#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
 
args <- commandArgs(trailingOnly = TRUE)

# Initialize a database with column names
chrom_summaries <- read.table(args[1], header = T, stringsAsFactors = F)
chrom_summaries <- chrom_summaries[FALSE,]

for (arg in args) {
	chrom_summaries <- rbind(chrom_summaries,
                             read.table(arg, header = T, stringsAsFactors = F))
}


write.table(chrom_summaries, file = "Chromosome_summary.txt", sep = "\t", row.names = FALSE) 
