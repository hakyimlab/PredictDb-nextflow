#! /usr/bin/env Rscript

# Force single-threaded BLAS before any BLAS-linked package loads. mclapply
# forks worker processes; a multi-threaded BLAS (OpenBLAS/MKL) can leave
# forked children with a corrupted internal thread pool (the OS threads that
# owned its locks don't exist post-fork), causing segfaults. Belt-and-suspenders
# alongside setting these in the calling shell, since env vars must be in
# place before the BLAS library initializes.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

suppressMessages(library(optparse))

# create options
option_list <- list(
  make_option(c("--chrom"), type="character", default=NULL,
              help="The chromosome you are processing",
              metavar="numeric"),
  make_option(c("--snp_annotation"), type="character", default=NULL,
              help="The snp annotation file for the chromosome you are processing",
              metavar="character"),
  make_option(c("--gene_annotation"), type="character", default=NULL,
              help="The gene annotation file for the chromosome you are processing",
              metavar="character"),
  make_option(c("--genotype_file"), type="character", default=NULL,
              help="The genotype file for the chromosome you are processing",
              metavar="character"),
  make_option(c("--gene_expression"), type="character", default=NULL,
              help="The gene expression file",
              metavar="character"),
  make_option(c("--covariates_file"), type="character", default=NULL,
              help="The covariates file",
              metavar="character"),
  make_option(c("--prefix"), type="character", default="elastic_net_cv",
              help="The prefix for the output files",
              metavar="character"),
  make_option(c("--nfolds"), type="numeric", default=5,
              help="The number of folds for the nested cross validation",
              metavar="numeric"),
  make_option(c("--n_times"), type="numeric", default=3,
              help="The number of repeated CV fold assignments to average over when selecting lambda",
              metavar="numeric"),
  make_option(c("--seed"), type="numeric", default=NA,
              help="Seed for reproducibility",
              metavar="numeric"),
  make_option(c("--n_cores"), type="numeric", default=1,
              help="Number of CPU cores to fit genes in parallel with (fork-based, Unix only). Intended to be set to the Nextflow process's task.cpus",
              metavar="numeric"),
  make_option(c("--chunk_size"), type="numeric", default=20,
              help="Number of genes to fit and write per batch, so results (esp. covariance data) for the whole chromosome are never all held in memory at once",
              metavar="numeric"))

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)


chrom <- args$chrom
snp_annot_file <- args$snp_annotation
gene_annot_file <- args$gene_annotation
genotype_file <- args$genotype_file
expression_file <- args$gene_expression
covariates_file <- args$covariates_file
prefix <- args$prefix
n_k_folds <- as.numeric(args$nfolds)
n_times <- as.numeric(args$n_times)
n_cores <- as.numeric(args$n_cores)
chunk_size <- as.numeric(args$chunk_size)

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(parallel))
"%&%" <- function(a,b) paste(a,b, sep = "")

# ---------------------------------------------------------------------------
# Data-loading / per-gene helper functions below are unchanged from
# elasticnet.R. Only main() is restructured to separate the random fold-id
# assignment (must stay sequential and in gene order to reproduce
# elasticnet.R's output) from the expensive cv.glmnet fitting (parallelized,
# since it is deterministic given its inputs and consumes no randomness).
# ---------------------------------------------------------------------------

get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F)

  if("refAllele" %in% names(snp_annot)) {snp_annot <- snp_annot %>% rename(ref_vcf = refAllele)}
  if("effectAllele" %in% names(snp_annot)) {snp_annot <- snp_annot %>% rename(alt_vcf = effectAllele)}

  snp_annot <- snp_annot %>%
    filter(!((ref_vcf == 'A' & alt_vcf == 'T') |
               (ref_vcf == 'T' & alt_vcf == 'A') |
               (ref_vcf == 'C' & alt_vcf == 'G') |
               (ref_vcf == 'G' & alt_vcf == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}


get_maf_filtered_genotype <- function(genotype_file_name,  maf, samples) {
  gt_df <- read.table(genotype_file_name, header = T, stringsAsFactors = F, row.names = 1)
  gt_df <- gt_df[,(colnames(gt_df) %in% samples )] %>% t() %>% as.data.frame()
  effect_allele_freqs <- colMeans(gt_df) / 2
  gt_df <- gt_df[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]
  gt_df
}

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding', 'pseudogene', 'lincRNA')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}

get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% t() %>% as.data.frame()
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  # Base-R subsetting instead of dplyr filter()/select(): with gt_df having
  # 100k+ columns, dplyr's tidyselect machinery re-scans/re-hashes all column
  # names on every call, which dominates runtime when this is called once per
  # gene. Base R's `[` uses a simple name-matched lookup instead. Row/column
  # order is preserved identically to the dplyr version, so output is unchanged.
  mask <- (snp_annot$pos >= (coords[1] - cis_window)) & (snp_annot$pos <= (coords[2] + cis_window)) & !is.na(snp_annot$rsid)
  snp_info <- snp_annot[mask, , drop = FALSE]
  if (nrow(snp_info) == 0)
    return(NA)
  matched_cols <- intersect(snp_info$varID, colnames(gt_df))
  if (length(matched_cols) == 0)
    return(NA) # none of the candidate varIDs exist in the gt_df dataset
  cis_gt <- gt_df[, matched_cols, drop = FALSE]
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  cov_df <- cov_df[,(names(cov_df) %in% samples )] %>% t() %>% as.data.frame()
  cov_df
}

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

do_elastic_net <- function(cis_gt, expr_adj, n_folds, cv_fold_ids, n_times, alpha) {
  cis_gt <- as.matrix(cis_gt)
  fit <- cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
  lambda_seq <- fit$lambda
  cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
  fits <- list()
  fits[[1]] <- fit
  cvms <- matrix(nrow = 100, ncol = n_times)
  cvms[1:length(fit$cvm),1] <- fit$cvm
  for (i in seq_len(n_times)[-1]) {
    fit <- cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
    fits[[i]] <- fit
    cvms[1:length(fit$cvm),i] <- fit$cvm
  }
  avg_cvm <- rowMeans(cvms)
  best_lam_ind <- which.min(avg_cvm)
  best_lambda <- lambda_seq[best_lam_ind]
  out <- list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda)
  out
}

evaluate_performance <- function(cis_gt, expr_adj, fit, best_lam_ind, best_lambda, cv_fold_ids, n_folds) {
  n_nonzero <- fit$nzero[best_lam_ind]
  if (n_nonzero > 0) {
    R2 <- rep(0, n_folds)
    corr_folds <- rep(0, n_folds)
    pval_folds <- rep(0, n_folds)
    zscore_folds <- rep(0, n_folds)

    for (j in (1:n_folds)) {
      fold_idxs <- which(cv_fold_ids[,1] == j)
      fitted_y <- fit$fit.preval[, best_lam_ind]

      tss <- sum(expr_adj[fold_idxs]**2)
      rss <- sum((expr_adj[fold_idxs] - fitted_y[fold_idxs])**2)
      R2[j] <- 1 - (rss/tss)

      corr_folds[j] <- ifelse(sd(fitted_y[fold_idxs]) != 0,
                       cor(expr_adj[fold_idxs], fitted_y[fold_idxs]), 0)
      pval_folds[j] <- ifelse(sd(fitted_y[fold_idxs]) != 0 & length(expr_adj[fold_idxs]) > 3,
                        cor.test(expr_adj[fold_idxs],fitted_y[fold_idxs])$p.value, 1)

      zscore_folds[j] <- atanh(corr_folds[j]) * sqrt(length(expr_adj[fold_idxs]) - 3) # Fisher transformation
    }

    best_fit <- fit$glmnet.fit
    expr_adj_pred <- predict(best_fit, as.matrix(cis_gt), s = best_lambda)
    tss <- sum(expr_adj**2)
    rss <- sum((expr_adj - expr_adj_pred)**2)

    n_samp <- length(expr_adj)
    weights <- best_fit$beta[which(best_fit$beta[,best_lam_ind] != 0), best_lam_ind]
    weighted_snps <- names(best_fit$beta[,best_lam_ind])[which(best_fit$beta[,best_lam_ind] != 0)]
    R2_mean <- mean(R2)
    R2_sd <- sd(R2)
    inR2 <- 1 - (rss/tss)
    rho_avg <- mean(corr_folds)
    rho_se <- sd(corr_folds)
    rho_avg_squared <- rho_avg**2

    # Stouffer's method for combining z scores.
    zscore_est <- sum(zscore_folds) / sqrt(n_folds)
    zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
    # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
    pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_folds, lower.tail = F)
    # Old way
    pred_perf <- summary(lm(expr_adj ~ fit$fit.preval[,best_lam_ind]))
    pred_perf_rsq <- pred_perf$r.squared

    out <- list(weights = weights, n_weights = n_nonzero, weighted_snps = weighted_snps, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se,
                rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, weighted_snps = NA, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pval_est = NA, rho_avg = NA, rho_se= NA,
                rho_zscore= NA, rho_avg_squared= NA, zscore_pval= NA)
  }
  out
}

compute_covariance <- function(gene_id, cis_gt, rsids, varIDs) {
  model_gt <- cis_gt[,varIDs, drop=FALSE]
  geno_cov <- cov(model_gt)
  n <- length(rsids)
  # Vectorized replacement for the old nested-loop + rbind()-per-pair pattern:
  # rbind() in a loop reallocates and copies the *entire* accumulated table on
  # every single call, which for a gene with many selected SNPs (n large, up
  # to O(n^2) pairs) can spike memory usage catastrophically and crash R
  # outright - this is what was segfaulting/aborting on genes with large
  # models. Builds the same (i, j) pairs, i=1..n, j=i..n, in one allocation.
  i_idx <- rep(seq_len(n), times = rev(seq_len(n)))
  j_idx <- unlist(lapply(seq_len(n), function(i) i:n))
  data.frame(gene = gene_id, rsid1 = rsids[i_idx], rsid2 = rsids[j_idx],
             covariance = geno_cov[cbind(i_idx, j_idx)])
}

# ---------------------------------------------------------------------------
# Phase A: sequential. Walks genes in the original order, building cis_gt /
# adjusted expression / cv_fold_ids exactly as elasticnet.R's main loop does.
# This is the *only* phase that touches the RNG (via sample() for fold ids),
# so it must run single-threaded and in-order to reproduce elasticnet.R's
# fold assignments bit-for-bit. It does not call cv.glmnet, so it is cheap.
# ---------------------------------------------------------------------------
prepare_gene_jobs <- function(genes, gene_annot, expr_df, gt_df, snp_annot, covariates_df,
                               cis_window, n_k_folds, n_times) {
  jobs <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)

    if (all(is.na(cis_gt))) {
      jobs[[i]] <- list(gene = gene, gene_name = gene_name, gene_type = gene_type, skip = TRUE)
      next
    }

    if (ncol(cis_gt) < 2) {
      # Matches elasticnet.R: n_snps_in_window is only ever populated when
      # ncol(cis_gt) >= 2, so it stays NA in the output here too.
      jobs[[i]] <- list(gene = gene, gene_name = gene_name, gene_type = gene_type, skip = TRUE)
      next
    }

    expression_vec <- expr_df[,i]
    if (is.null(covariates_df)) {
      adj_expression <- as.matrix(expression_vec)
      rownames(adj_expression) <- rownames(expr_df)
    } else {
      adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
    }
    # Filter out samples that are in both expression and genotype data.
    adj_expression <- as.matrix(adj_expression[(rownames(adj_expression) %in% rownames(cis_gt)),])
    # sort row names to be in the same order for X and y
    cis_gt <- cis_gt[match(rownames(adj_expression), rownames(cis_gt)),]

    # make the cv folds (the only randomness in the whole pipeline)
    n_filt <- length(adj_expression)
    cv_fold_ids <- matrix(nrow = n_filt, ncol = n_times)
    for (j in 1:n_times)
      cv_fold_ids[,j] <- sample(1:n_k_folds, n_filt, replace = TRUE)

    jobs[[i]] <- list(gene = gene, gene_name = gene_name, gene_type = gene_type, skip = FALSE,
                       coords = coords, adj_expression = adj_expression, cv_fold_ids = cv_fold_ids,
                       n_snps_in_window = ncol(cis_gt))
  }
  jobs
}

# ---------------------------------------------------------------------------
# Phase B: parallel. Recomputes cis_gt (cheap column slice) from the shared
# gt_df/snp_annot every fork already has via copy-on-write, then fits the
# model. Deterministic given its inputs - consumes no randomness - so it is
# safe to run in any order across any number of workers.
# ---------------------------------------------------------------------------
run_gene_job <- function(job, gt_df, snp_annot, snp_annot_lookup, cis_window, n_k_folds, n_times, alpha) {
  if (job$skip) return(job)

  cis_gt <- get_cis_genotype(gt_df, snp_annot, job$coords, cis_window)
  cis_gt <- cis_gt[match(rownames(job$adj_expression), rownames(cis_gt)),]

  out <- tryCatch({
    elnet_out <- do_elastic_net(cis_gt, job$adj_expression, n_k_folds, job$cv_fold_ids, n_times, alpha)
    eval <- evaluate_performance(cis_gt, job$adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind,
                                  elnet_out$best_lambda, job$cv_fold_ids, n_k_folds)

    weights_info <- NULL
    covariance_df <- NULL
    if (!is.na(eval$n_weights) && eval$n_weights > 0) {
      weights_info <- snp_annot_lookup %>% filter(varID %in% eval$weighted_snps) %>% select(rsid, varID, ref_vcf, alt_vcf)
      weights_info$gene <- job$gene
      weights_info <- weights_info %>% merge(data.frame(weights = eval$weights, varID = eval$weighted_snps), by = 'varID') %>%
                        select(gene, rsid, varID, ref_vcf, alt_vcf, weights)
      covariance_df <- compute_covariance(job$gene, cis_gt, weights_info$rsid, weights_info$varID)
    }

    c(job, list(elnet_out = elnet_out, eval = eval, weights_info = weights_info,
                covariance_df = covariance_df, error = NULL))
  }, error = function(cond) {
    c(job, list(elnet_out = NULL, eval = NULL, weights_info = NULL, covariance_df = NULL,
                error = paste('gene', job$gene, ':', conditionMessage(cond))))
  })
  out
}

write_gene_result <- function(res, alpha, model_summary_file, weights_file, covariance_file) {
  model_summary <- c(res$gene, res$gene_name, res$gene_type, alpha, NA, NA, NA, NA, 0,
                      NA, NA, NA, NA, NA, NA, NA, NA, NA)

  if (!is.null(res$error)) {
    cat("Warning:", res$error, "\n")
  } else if (!res$skip && !is.null(res$eval)) {
    model_summary <- c(res$gene, res$gene_name, res$gene_type, alpha, res$n_snps_in_window,
                       res$elnet_out$min_avg_cvm, res$elnet_out$best_lam_ind, res$elnet_out$best_lambda,
                       res$eval$n_weights, res$eval$R2_mean, res$eval$R2_sd, res$eval$inR2,
                       res$eval$pval_est, res$eval$rho_avg, res$eval$rho_se, res$eval$rho_zscore,
                       res$eval$rho_avg_squared, res$eval$zscore_pval)
    if (!is.null(res$weights_info)) {
      write.table(res$weights_info, file = weights_file, append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE, sep = '\t')
      write.table(res$covariance_df, file = covariance_file, append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE, sep = ' ')
    }
  }
  # else: res$skip == TRUE (no SNPs in window, or <2 SNPs) -> keep the default NA row,
  # matching elasticnet.R exactly.
  write(model_summary, file = model_summary_file, append = TRUE, ncol = 18, sep = '\t')
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,
                 covariates_file, chrom, prefix, maf=0.01, n_times=3, n_k_folds=10,
                 seed=NA, cis_window=1e6, alpha=0.5, n_cores=1, chunk_size=200) {
  # Read in data----
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df <- get_gene_expression(expression_file, gene_annot)
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(expr_df)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)

  covariates_df <- NULL
  if (!is.null(covariates_file)) {
    covariates_df <- get_covariates(covariates_file, samples)
  }

  # Set seed and cross-validation fold ids----
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)

  # Prepare output data----
  model_summary_file <- './summary/' %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
  model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'cv_mse',
                          'lambda_iteration', 'lambda_min', 'n_snps_in_model',
                          'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2', 'cv_fisher_pval',
                          'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval')
  write(model_summary_cols, file = model_summary_file, ncol = 18, sep = '\t')

  weights_file <- './weights/' %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'
  weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
  write(weights_col, file = weights_file, ncol = 6, sep = '\t')

  tiss_chr_summ_f <- './chrom_summary/' %&% prefix %&% '_chr' %&% chrom %&% '_summary.txt'
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')

  covariance_file <- './covariances/' %&% prefix %&% '_chr' %&% chrom %&% '_covariances.txt'
  covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')

  # Phase A: sequential fold-id assignment (the only randomness in the pipeline)
  cat("Preparing", n_genes, "genes (sequential)...\n")
  jobs <- prepare_gene_jobs(genes, gene_annot, expr_df, gt_df, snp_annot, covariates_df,
                            cis_window, n_k_folds, n_times)

  # Phase B+C: parallel model fitting, chunked so that at most chunk_size genes'
  # worth of results (weights + covariance data frames, which can be large for
  # genes with many SNPs in their model) are ever held in memory at once. Each
  # chunk is fit in parallel, written to disk immediately in original gene
  # order, then discarded before moving to the next chunk.
  n_chunks <- ceiling(n_genes / chunk_size)
  cat("Fitting", n_genes, "genes across", n_cores, "core(s) in", n_chunks, "chunk(s) of up to", chunk_size, "genes...\n")
  for (chunk_idx in seq_len(n_chunks)) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_genes)
    cat("Chunk", chunk_idx, "/", n_chunks, "(genes", start_idx, "-", end_idx, ")\n")

    # mc.preschedule=FALSE: dispatch one gene at a time to whichever worker is
    # free next, instead of splitting the chunk into fixed blocks up front.
    # Per-gene cost varies enormously with local SNP density (a few hundred to
    # several thousand cis-SNPs), so static blocks can strand one worker with
    # several expensive genes while others sit idle.
    chunk_results <- mclapply(jobs[start_idx:end_idx], run_gene_job, gt_df = gt_df, snp_annot = snp_annot,
                               snp_annot_lookup = snp_annot, cis_window = cis_window,
                               n_k_folds = n_k_folds, n_times = n_times, alpha = alpha,
                               mc.cores = n_cores, mc.preschedule = FALSE)

    for (res in chunk_results) {
      write_gene_result(res, alpha, model_summary_file, weights_file, covariance_file)
    }

    rm(chunk_results)
    gc()
  }
}

#Run analysis
main(snp_annot_file, gene_annot_file, genotype_file, expression_file,
    covariates_file, as.numeric(chrom), prefix, n_times = n_times, n_k_folds = n_k_folds,
    seed = args$seed, n_cores = n_cores, chunk_size = chunk_size)
