#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
 
args <- commandArgs(trailingOnly = TRUE)

# Initialize a database with column names
cova_df <- read.table(args[1], header = T, stringsAsFactors = F)
cova_df <- cova_df[FALSE,]

for (arg in args) {
	cova_df <- rbind(cova_df,
                             read.table(arg, header = T, stringsAsFactors = F))
}

write.table(cova_df, file = "Covariances.txt", sep = "\t", row.names = FALSE, quote = FALSE) 
