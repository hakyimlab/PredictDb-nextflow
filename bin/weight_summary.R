#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
 
args <- commandArgs(trailingOnly = TRUE)

# Initialize a database with column names
weights <- read.table(args[1], header = T, stringsAsFactors = F)
weights <- weights[FALSE,]

for (arg in args) {
	weights <- rbind(weights,
                             read.table(arg, header = T, stringsAsFactors = F))
}


write.table(weights, file = "Weight_summary.txt", sep = "\t", row.names = FALSE) 
