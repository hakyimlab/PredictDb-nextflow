#! /usr/bin/env nextflow

/*
==================================================
                    PredictDb
==================================================
 PredictDb Analysis Pipeline.
 #### Homepage / Documentation
 https:://github.com/hakyimlab/PredictDb-nextflow
---------------------------------------------------
*/


//Enable DSL 2 syntax
nextflow.enable.dsl = 2


def helpMessage() {
    log.info """
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --gene_annotation 'gene_annot.gtf' --snp_annotation 'snp_annnotation_file.vcf' --genotype 'genotype_file.txt' --gene_exp 'Normalized_gene expression.txt'
    

    Mandatory arguments:
       --gene_annotation [file]		Path to the gene annotation file .gtf format (must be surrounded with quotes)
       --snp_annotation [file]		Path to the SNP annotation file .vcf or .txt format (must be surrounded with quotes)
       --genotype [file]			Path to the dosage file .txt format and TAB separated (must be surrounded with quotes)
       --gene_exp [file]			Path to the gene expression file .csv,.txt or .tsv format (must be surrounded with quotes)

    Options:
       --pca [bool]				Specifies if you want to perform PCA on the gene expression data. The default is false
       --peer [bool]                        Specifies if you want to perform PEER on the gene expression data. The default is false
       --covariates [file]			Specifies if you have covariates to be regressed out on the gene expression data.
       --nfolds [num]			The number of folds to split your data into. The default is 10
       --nested_cv [bool]			Run the nested cross validated elasticnet, default: cross validated elasticnet.
       --keepIntermediate [bool]		Specifies if you want to keep intermediate files
       --prefix [str]			The prefix of your output files, we recommend using the population name used in the training. If not provided default name is used
       --outdir [file]			The output directory where the results will be saved
       --email [email]			Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
       --email_on_fail [email]		Same as --email, except only send mail if the workflow is not successful
       -name [str]				Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// catching both -name and --name if specified by user
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

//Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']              = custom_runName ?: workflow.runName
summary['Gene annotation']       = params.gene_annotation
summary['SNP annotation']        = params.snp_annotation
summary['Expression file']       = params.gene_exp
summary['Genotype file']         = params.genotype
if (params.covariates) {
    summary['Covariates']    = params.covariates
}
summary['Max Resources']         = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']            = params.outdir
summary['Launch dir']            = workflow.launchDir
summary['Working dir']           = workflow.workDir
summary['Script dir']            = workflow.projectDir
summary['User']                  = workflow.userName
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



// Set files path
gtf = Channel
    .fromPath(params.gene_annotation, checkIfExists: true)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gene_annotation}" }

snp_annot = Channel
    .fromPath(params.snp_annotation, checkIfExists: true)
    .ifEmpty { exit 1, "SNP annotation file not found: ${params.snp_annotation}" }

genotype = Channel
    .fromPath(params.genotype, checkIfExists: true)
    .ifEmpty { exit 1, "Genotype file not found: ${params.genotype}" }

geneExp = Channel
    .fromPath(params.gene_exp, checkIfExists: true)
    .ifEmpty { exit 1, "Gene expression file not found: ${params.gene_exp}" }

p_covs = null
if (params.covariates) {
    p_covs = Channel
        .fromPath(params.covariates, checkIfExists: true)
        .ifEmpty { exit 1, "Covariates file not found: ${params.covariates}" }
}


// import workflow
include { UNZIP;
            PREPROCESS;
            COVS;
            COMBINE_COVS;
            TRAIN_MODEL;
            CREATE_DB } from './workflow/main_workflow.nf' addParams(outdir: "${params.outdir}")


workflow {
    main:
      UNZIP(gtf,snp_annot,genotype,geneExp)
      PREPROCESS(UNZIP.out.gtf,UNZIP.out.snp_annot,UNZIP.out.genotype,UNZIP.out.geneExp)
      COVS(PREPROCESS.out.tr_expr_peer,PREPROCESS.out.tr_expr_count)
      COMBINE_COVS(COVS.out.covs,p_covs)
      TRAIN_MODEL(PREPROCESS.out.gene_annot,PREPROCESS.out.snp_files,
                  PREPROCESS.out.genotype_files,PREPROCESS.out.tr_expr,COMBINE_COVS.out.covs)
      CREATE_DB(TRAIN_MODEL.out.chrom_summaries.collect(),
                TRAIN_MODEL.out.weight_summaries.collect(),
                TRAIN_MODEL.out.model_summaries.collect(),
                TRAIN_MODEL.out.model_covariances.collect())
}

