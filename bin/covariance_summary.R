#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_file <- "Covariances.txt"

open_maybe_gz <- function(path, mode = "r") {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, mode) else file(path, mode)
}

# Process one chromosome's covariance file at a time, in full, instead of
# loading every chromosome into data.frames and rbind()-ing them all together
# (which held the whole genome-wide covariance table in memory at once - the
# same rbind()-growth problem already fixed in fast_elasticnet.R's
# compute_covariance, but here compounded across every chromosome file
# simultaneously - and crashed with a bus error once there were many SNPs).
# A single chromosome's file is bounded in size, so reading it whole is safe
# given this step's memory budget, and avoids the overhead of splitting each
# file into many smaller read/write calls. Only the header from the first
# file is kept, and every other file's header is checked against it rather
# than relied upon implicitly (the way rbind()'s column matching used to).
out_con <- file(out_file, open = "w")
first_header <- NULL

for (i in seq_along(args)) {
  in_con <- open_maybe_gz(args[i])
  all_lines <- readLines(in_con)
  close(in_con)

  header <- gsub("[ \t]+", "\t", all_lines[1])
  if (i == 1) {
    first_header <- header
    writeLines(header, out_con)
  } else if (!identical(header, first_header)) {
    stop("Column header in ", args[i], " does not match ", args[1])
  }

  if (length(all_lines) > 1) {
    writeLines(gsub("[ \t]+", "\t", all_lines[-1]), out_con)
  }
  rm(all_lines)
}

close(out_con)
