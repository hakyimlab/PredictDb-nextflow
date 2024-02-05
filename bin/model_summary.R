#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
 
args <- commandArgs(trailingOnly = TRUE)

# Initialize a database with column names
model_summaries <- read.table(args[1], header = T, stringsAsFactors = F, fill = T)
model_summaries <- model_summaries[FALSE,]

for (arg in args) {
	model_summaries <- rbind(model_summaries,
                             read.table(arg, header = T, stringsAsFactors = F,fill =T))
}


write.table(model_summaries, file = "Model_summary.txt", sep = "\t", row.names = FALSE) 
