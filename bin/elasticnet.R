#! /usr/bin/env Rscript
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
  make_option(c("--seed"), type="numeric", default=1824,
              help="Seed for reproducibility",
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

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(tools))
"%&%" <- function(a,b) paste(a,b, sep = "")


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
    filter((chr == chrom) & gene_type %in% gene_types) %>% 
    mutate(gene_id = file_path_sans_ext(gene_id))  # remove gene ver from gene names
  gene_df
}

get_gene_type <- function(gene_annot, gene) {
  gene_annot <- gene_annot %>% mutate(gene_id = file_path_sans_ext(gene_id))
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% t() %>% as.data.frame()
  colnames(expr_df) <- file_path_sans_ext(colnames(expr_df)) # remove gene version number
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  #Check of the varID exist in the data
  if (TRUE %in% (snp_info$varID %in% names(gt_df))) {
    cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
  } else {
    return(NA) # the varID doesn't exist in the gt_df dataset
  }
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
  #tryCatch({
    cis_gt <- as.matrix(cis_gt)
    fit <- cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
    lambda_seq <- fit$lambda
    cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
    fits <- list()
    fits[[1]] <- fit
    cvms <- matrix(nrow = 100, ncol = n_times)
    cvms[1:length(fit$cvm),1] <- fit$cvm
    for (i in 2:(n_times)) {
      fit <- cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
      fits[[i]] <- fit
      cvms[1:length(fit$cvm),i] <- fit$cvm
    }
    avg_cvm <- rowMeans(cvms)
    best_lam_ind <- which.min(avg_cvm)
    best_lambda <- lambda_seq[best_lam_ind]
    out <- list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda)
  #  out
  #},error = function(cond) {
  #  message('Error')
  #  message(geterrmessage())
  #  out <- list()
  #  out
  #}
  #)
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
    #f_stat <- ((tss - rss) / n_nonzero) / (rss / (n_samp - ncol(cis_gt) - 1))
    #p_val <- pf(f_stat, n_samp - 1, n_samp - ncol(cis_gt) - 1, lower.tail = FALSE)
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
    
    #one_sided_pval <- cor.test(expr_adj, fit$fit.preval[,best_lam_ind], alternative = 'greater')$p.value
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

do_covariance <- function(gene_id, cis_gt, rsids, varIDs, out_file) {
  model_gt <- cis_gt[,varIDs, drop=FALSE]
  geno_cov <- cov(model_gt)
  #print(dim(geno_cov))
  cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
  for (i in 1:length(rsids)) {
    for (j in i:length(rsids)) {
      #print(c(i, j))
      cov_df <- tryCatch(rbind(cov_df, data.frame(gene=gene_id,rsid1=rsids[i], rsid2=rsids[j], covariance=geno_cov[i,j])),
                         error = function(cond) browser())
    }
  }
  write.table(cov_df, file = out_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,
                 covariates_file, chrom, prefix, maf=0.01, n_times=3, n_k_folds=10,
                 seed=NA, cis_window=1e6, alpha=0.5) {
  # Read in data----
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df <- get_gene_expression(expression_file, gene_annot)
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(expr_df)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)
  
  if (!is.null(args$covariates_file)) {
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
  
  # Attempt to build model for each gene----
  for (i in 1:n_genes) {
    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    model_summary <- c(gene, gene_name, gene_type, alpha, NA, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)

    if (all(is.na(cis_gt))) {
      # No snps within window for gene.
      write(model_summary, file = model_summary_file, append = TRUE, ncol = 18, sep = '\t')
      next
    }

    if (ncol(cis_gt) >= 2) {
      expression_vec <- expr_df[,i]
      if (is.null(args$covariates_file)) {
        adj_expression <- as.matrix(expression_vec)
        rownames(adj_expression) <- rownames(expr_df)
      } else {
        adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
      }
      # Filter out samples that are in both expression and genotype data.
      adj_expression <- as.matrix(adj_expression[(rownames(adj_expression) %in% rownames(cis_gt)),])
      # sort row names to be in the same order for X and y
      cis_gt = cis_gt[match(rownames(adj_expression), rownames(cis_gt)),]

      # make the cv folds
      n_filt <- length(adj_expression)
      cv_fold_ids <- matrix(nrow = n_filt, ncol = n_times)
      for (j in 1:n_times)
        cv_fold_ids[,j] <- sample(1:n_k_folds, n_filt, replace = TRUE)

      elnet_out <- do_elastic_net(cis_gt, adj_expression, n_k_folds, cv_fold_ids, n_times, alpha)
      if (length(elnet_out) > 0) {
        eval <- evaluate_performance(cis_gt, adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind, elnet_out$best_lambda, cv_fold_ids, n_k_folds)
        model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), elnet_out$min_avg_cvm, elnet_out$best_lam_ind,
                           elnet_out$best_lambda, eval$n_weights, eval$R2_mean, eval$R2_sd, eval$inR2,
                           eval$pval_est, eval$rho_avg, eval$rho_se, eval$rho_zscore, eval$rho_avg_squared, eval$zscore_pval)
          #return(list(eval = eval, model_summary = model_summary, elnet_out = elnet_out, snp_annot = snp_annot))
        if (eval$n_weights > 0) {
          weighted_snps_info <- snp_annot %>% filter(varID %in% eval$weighted_snps) %>% select(rsid, varID, ref_vcf, alt_vcf)
          if (nrow(weighted_snps_info) == 0)
            browser()
          weighted_snps_info$gene <- gene
          weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights = eval$weights, varID=eval$weighted_snps), by = 'varID') %>%
                                  select(gene, rsid, varID, ref_vcf, alt_vcf, weights)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          do_covariance(gene, cis_gt, weighted_snps_info$rsid, weighted_snps_info$varID, covariance_file)
        }
      }
    }
    write(model_summary, file = model_summary_file, append = TRUE, ncol = 18, sep = '\t')
  }
}

#Run analysis
main(snp_annot_file, gene_annot_file, genotype_file, expression_file, 
    covariates_file, as.numeric(chrom), prefix,n_k_folds = n_k_folds)
