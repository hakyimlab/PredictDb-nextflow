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
        --gtf
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

// Validate input files
if (params.gtf) {
 if (hasExtension(params.gtf, 'gz')) {
  gtf_gz = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  } else {
  newch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  }
}


// Gunzip the gtf file
if (params.gtf && hasExtension(params.gtf, 'gz')) {
    process gunzip_gtf {
        tag "$gz"
        publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
                   saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

        input:
        file gz from gtf_gz

        output:
        file "${gz.baseName}" into newch

        script:
        """
        gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
        """
    }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
