#! /usr/bin/env nextflow

/*
========================================================================================
                         predixcan
========================================================================================
 PrediXcan Analysis Pipeline.
 #### Homepage / Documentation
 https:
----------------------------------------------------------------------------------------
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
    path "gene_annot.parsed.txt" into ch

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
    path "snp_annot.*.txt" into ch2

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
    path "genotype.*.txt" into ch3

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
    path "transposed_gene_exp.csv" into ch4

    script:
    """
    transpose_gene_expr.R ${geneExp} transposed_gene_exp.csv
    """
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
