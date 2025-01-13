#! /usr/bin/env nextflow

process model_training_w_covs {
    tag "model training"
    label "process_high"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/models" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    
    input:
    file covariates 
    file expression 
    file gene_annot 
    tuple val(chrom), file(snp_file2:'snp_file'), file('genotype_file') 

    output:
    path "weights/*", emit: weight_summaries
    path "summary/*", emit: model_summaries
    path "chrom_summary/*", emit: chrom_summaries
    path "covariances/*", emit: all_covariances

    script:
    prefix = params.prefix
    nfolds = params.nfolds
    """
    mkdir -p summary weights covariances chrom_summary
    elasticnet.R \
        --chrom $chrom \
        --snp_annotation snp_file \
        --gene_annotation $gene_annot \
        --genotype_file genotype_file \
        --gene_expression $expression \
        --covariates_file $covariates \
        --prefix $prefix \
        --nfolds $nfolds
    """
}


process model_training_wo_covs {
    tag "model training"
    label "process_high"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/models" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    
    input:
    file expression 
    file gene_annot
    tuple val(chrom), file(snp_file2:'snp_file'), file('genotype_file')

    output:
    path "weights/*", emit: weight_summaries
    path "summary/*", emit: model_summaries
    path "chrom_summary/*", emit: chrom_summaries
    path "covariances/*", emit: all_covariances

    script:
    prefix = params.prefix
    nfolds = params.nfolds
    """
    mkdir -p summary weights covariances chrom_summary
    elasticnet.R \
        --chrom $chrom \
        --snp_annotation snp_file \
        --gene_annotation $gene_annot \
        --genotype_file genotype_file \
        --gene_expression $expression \
        --prefix $prefix \
        --nfolds $nfolds
    """
}

// nected cv 
process model_training_w_covs_nested {
    tag "model training"
    label "process_high"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/models" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    
    input:
    file covariates 
    file expression 
    file gene_annot 
    tuple val(chrom), file(snp_file2:'snp_file'), file('genotype_file') 

    output:
    path "weights/*", emit: weight_summaries
    path "summary/*", emit: model_summaries
    path "chrom_summary/*", emit: chrom_summaries
    path "covariances/*", emit: all_covariances

    script:
    prefix = params.prefix
    nfolds = params.nfolds
    """
    mkdir -p summary weights covariances chrom_summary
    nested_cv_elasticnet.R \
        --chrom $chrom \
        --snp_annotation snp_file \
        --gene_annotation $gene_annot \
        --genotype_file genotype_file \
        --gene_expression $expression \
        --covariates_file $covariates \
        --prefix $prefix \
        --nfolds $nfolds
    """
}


process model_training_wo_covs_nested {
    tag "model training"
    label "process_high"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/models" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    
    input:
    file expression 
    file gene_annot
    tuple val(chrom), file(snp_file2:'snp_file'), file('genotype_file')

    output:
    path "weights/*", emit: weight_summaries
    path "summary/*", emit: model_summaries
    path "chrom_summary/*", emit: chrom_summaries
    path "covariances/*", emit: all_covariances

    script:
    prefix = params.prefix
    nfolds = params.nfolds
    """
    mkdir -p summary weights covariances chrom_summary
    nested_cv_elasticnet.R \
        --chrom $chrom \
        --snp_annotation snp_file \
        --gene_annotation $gene_annot \
        --genotype_file genotype_file \
        --gene_expression $expression \
        --prefix $prefix \
        --nfolds $nfolds
    """
}
