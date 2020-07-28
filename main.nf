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

    nextflow run main.nf --vcf 'vcffile.vcf' --gene_annot 'gene_annnotation_file'
    

    Mandatory arguments:
        --vcf
        --gene_annot
        --snp_annot
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

