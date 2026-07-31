#! /usr/bin/env Rscript

suppressMessages(library(RSQLite))

args <- commandArgs(trailingOnly = TRUE)

filtered_db <- args[1]
unfiltered_cov <- args[2]
filtered_cov <- args[3]

open_maybe_gz <- function(path, mode = "r") {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, mode) else file(path, mode)
}

# load db
driver <- dbDriver("SQLite")
in_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from extra')
dbDisconnect(in_conn)

keep_genes <- unique(model_summaries$gene)

# Stream the covariance file (the gzipped output of covariance_summary.R) in
# bounded-size batches instead of loading it whole into a data.frame - with
# many SNPs the combined genome-wide covariance file can be too large to hold
# in memory at once, the same failure mode fixed in covariance_summary.R.
# Only the first (GENE) field of each line needs to be inspected to decide
# whether to keep it. Measured peak working set (raw lines + extracted gene
# names + filtered output, held at once per batch) is ~230MB per 1e6 lines,
# so 2e7 lines/batch is ~4.6GB peak - safe with plenty of headroom for R's
# own overhead under either a 16GB or 32GB allocation for this step.
batch_size <- 2e7
in_con <- open_maybe_gz(unfiltered_cov)
out_con <- file(filtered_cov, open = "w")

header <- readLines(in_con, n = 1)
writeLines(header, out_con)

repeat {
  lines <- readLines(in_con, n = batch_size)
  if (length(lines) == 0) break
  genes <- sub("[ \t].*", "", lines)
  keep <- genes %in% keep_genes
  if (any(keep)) writeLines(lines[keep], out_con)
}

close(in_con)
close(out_con)
