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

def helpMessage() {
    log.info """
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --gtf 'vcffile.vcf' --gene_annot 'gene_annnotation_file'
    

    Mandatory arguments:
        --gtf
        --snp
        --genotype
        --gene_exp
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate input files
if (params.gtf) {
 if (hasExtension(params.gtf, 'gz')) {
  gtf_gz = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  } else {
  gene_annot = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  }
}

if (params.snp) {
 if (hasExtension(params.snp, 'gz')) {
  snp_gz = Channel
        .fromPath(params.snp, checkIfExists: true)
        .ifEmpty { exit 1, "SNP annotation file not found: ${params.snp}" }
  } else {
  snp_annot = Channel
        .fromPath(params.snp, checkIfExists: true)
        .ifEmpty { exit 1, "SNP annotation file not found: ${params.snp}" }
  }
}

if (params.genotype) {
 if (hasExtension(params.genotype, 'gz')) {
  genotype_gz = Channel
        .fromPath(params.genotype, checkIfExists: true)
        .ifEmpty { exit 1, "Genotype file not found: ${params.genotype}" }
  } else {
  gene_split = Channel
        .fromPath(params.genotype, checkIfExists: true)
        .ifEmpty { exit 1, "Genotype file not found: ${params.genotype}" }
  }
}

if (params.gene_exp) {
 if (hasExtension(params.gene_exp, 'gz')) {
  geneExp_gz = Channel
        .fromPath(params.gene_exp, checkIfExists: true)
        .ifEmpty { exit 1, "Gene expression file not found: ${params.gene_exp}" }
  } else {
  gene_expr = Channel
        .fromPath(params.gene_exp, checkIfExists: true)
        .ifEmpty { exit 1, "Gene expression file not found: ${params.gene_exp}" }
  }
}

/*
 * -------------------------------------------------
 *  Gunzip all the inputs if they are gzipped (.gz)
 * -------------------------------------------------
 */

if (params.gtf && hasExtension(params.gtf, 'gz')) {
    process gunzip_gtf {
        tag "$gz"
        publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
                   saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
        input:
        path gz from gtf_gz

        output:
        path "${gz.baseName}" into gene_annot

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}

if (params.snp && hasExtension(params.snp, 'gz')) {
    process gunzip_annot {
        tag "$gz"
        publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
                   saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
        input:
        path gz from snp_gz

        output:
        path "${gz.baseName}" into snp_annot

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}

if (params.genotype && hasExtension(params.genotype, 'gz')) {
    process gunzip_genotype {
        tag "$gz"
        publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
                   saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
        input:
        path gz from genotype_gz

        output:
        path "${gz.baseName}" into gene_split

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}

if (params.gene_exp && hasExtension(params.gene_exp, 'gz')) {
    process gunzip_geneExpression {
        tag "$gz"
        publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
                   saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
        input:
        path gz from geneExp_gz

        output:
        path "${gz.baseName}" into gene_expr

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}

/*
 * -------------------------------------------------
 *  Pre-process all the input files
 * -------------------------------------------------
 */

process gene_annotation {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/gtf" : false },    
	       saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path gtf from gene_annot

    output:
    path "gene_annot.parsed.txt" into parsed_annot

    script:
    """
    parse_gtf.py ${gtf} gene_annot.parsed.txt
    """
}

process split_annotation {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/snp_annot" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path snp from snp_annot

    output:
    path "snp_annot.*.txt" into snp_files

    script:
    """
    split_snp_annot_by_chr.py ${snp} snp_annot
    """
}

process split_genotype {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/genotype" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path genot from gene_split

    output:
    path "genotype.*.txt" into genotype_files

    script:
    """
    split_genotype_by_chr.py ${genot} genotype
    """
}

process transpose_geneExpression {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/transpose_gene-expr" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path geneExp from gene_expr

    output:
    path "transposed_gene_exp.csv" into tr_expr, gene_cols 

    script:
    """
    transpose_gene_expr.R ${geneExp} transposed_gene_exp.csv 
    wc -l transposed_gene_exp.csv | cut -d " " -f1 > count.txt
    """
}

/* Calculate PEER factors
 * ------------------------
 * If the number of samples is greater than or equal to 350, we use 60 PEER factors, 
 * If the number of samples is between 250 and 350, we use 45 PEER factors,
 * If the number of samples is between 150 and 250, we use 30 PEER factors, 
 * and if the number of samples is less than 150 we use 15 PEER factors.
 */

process generate_peer_factors {
    tag "PEER Factors"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/PEER" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path csv from tr_expr

    output:
    path "calculated_peers/X.csv" into peers

    script:
    """
    peertool -f ${csv} -n 5 --has_header -o calculated_peers
    """
}

process linear_regression {
    tag "regression"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/Covariates" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path peer from peers
    path gene_expr from gene_cols

    output:
    path "covariates.txt" into covariate_file
    path "transformed_expression.txt" into final_expr
    
    script:
    """
    process_covariates.R ${peer} ${gene_expr} covariates.txt transformed_expression.txt
    """
}

// Map files for each chromosome together
snp_files
     .flatMap()
     .map { file -> tuple(getChromID(file), file) }
     .set { map_snp }

genotype_files
     .flatMap()
     .map { file -> tuple(getChromID(file), file) }
     .set { map_genotype }

map_snp.join(map_genotype).set {snp_genotype_files}

//snp_genotype_files.subscribe onNext: { println it }


process model_training {
    tag "training"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/models" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    
    input:
    file covariates from covariate_file.first()
    file expression from final_expr.first()
    file gene_annot from parsed_annot.first()
    tuple val(chrom), file(snp_file2:'snp_file'), file('genotype_file') from snp_genotype_files

    output:
    path "weights/*" into weight_summaries
    path "summary/*" into model_summaries
    path "chrom_summary/*" into chrom_summaries
    path "covariances/*" into all_covariances

    script:
    prefix = params.prefix
    """
    mkdir -p summary weights covariances chrom_summary
    gtex_v7_nested_cv_elnet.R $chrom snp_file $gene_annot genotype_file $expression $covariates $prefix
    """
}

/*
ccchq.flatMap()
     .collect()
     .subscribe onNext: { println it }


ccchq
     .map { file -> tuple(getTrainID(file[0]), file[0], file[1]) }
     .toSortedList({ a, b -> a[0] <=> b[0] })
     .set { summaries }

*/
process collectModel_summaries {
    tag "database"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'

    input:
    path model from model_summaries.collect()

    output:
    path "Model_summary.txt" into all_model_sum

    script:
    """
    model_summary.R $model*
    """
}

process collectWeight_summaries {
    tag "database"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'

    input:
    path weight from weight_summaries.collect()

    output:
    path "Weight_summary.txt" into all_weight_sum

    script:
    """
    weight_summary.R $weight*
    """
}

process collectChrom_summaries {
    tag "database"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'

    input:
    path chrom from chrom_summaries.collect()

    output:
    path "Chromosome_summary.txt" into all_chrom_sum

    script:
    """
    chrom_summary.R $chrom*
    """
}

process make_database {
    tag "database"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'

    input:
    path models from all_model_sum
    path weights from all_weight_sum
    path chroms from all_chrom_sum

    output:

    script:
    """
    
    """
}


// Map the model training files together
def getTrainID( file ) {
    file.name.toString().find(/(.chr)(\d+)(_model_)/) { match, pref, chrom, ext -> chrom }
}

// Get the chromosome number
def getChromID( file ) {
    file.name.toString().find(/(.chr)(\d+)(.txt)/) { match, pref, chrom, ext -> chrom }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

