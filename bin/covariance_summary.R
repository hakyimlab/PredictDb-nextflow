#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_file <- "Covariances.txt"

open_maybe_gz <- function(path, mode = "r") {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, mode) else file(path, mode)
}

# Stream each chromosome's covariance file in bounded-size batches instead of
# reading it whole - some chromosomes (e.g. chr1) can he hige hence running to
# meory issues and crushing. Only the header from the first file is kept, and every other
# file's header is checked against it rather than relied upon implicitly (the
# way rbind()'s column matching used to).

batch_size <- 4e7

normalize_sep <- function(x) chartr(" ", "\t", x)

out_con <- file(out_file, open = "w")
first_header <- NULL

for (i in seq_along(args)) {
  in_con <- open_maybe_gz(args[i])

  header <- normalize_sep(readLines(in_con, n = 1))
  if (i == 1) {
    first_header <- header
    writeLines(header, out_con)
  } else if (!identical(header, first_header)) {
    stop("Column header in ", args[i], " does not match ", args[1])
  }

  repeat {
    lines <- readLines(in_con, n = batch_size)
    if (length(lines) == 0) break
    writeLines(normalize_sep(lines), out_con)
  }

  close(in_con)
}

close(out_con)
