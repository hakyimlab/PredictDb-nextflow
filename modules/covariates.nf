#! /usr/bin/env nextflow

process generate_peer_factors {
    tag "PEER Factors"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/PEER" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input: 
    path csv
    tcount

    output:
    path "calculated_peers/X.csv", emit: peers

    script:
    no_samples = tcount.toInteger() - 1
    if( no_samples > 0 && no_samples < 150 )
        """
        peertool -f ${csv} -n 15 --has_header -o calculated_peers
        """

    else if( no_samples >= 150 && no_samples < 250 )
        """
        peertool -f ${csv} -n 30 --has_header -o calculated_peers
        """

    else if( no_samples >= 250 && no_samples < 350 )
        """
        peertool -f ${csv} -n 45 --has_header -o calculated_peers
        """

    else if( no_samples >= 350 )
        """
        peertool -f ${csv} -n 60 --has_header -o calculated_peers
        """

    else
        error "Invalid number of samples in gene expression file: ${no_samples}"
}

process process_peer_factors {
    tag "PEER Factors"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/Covariates" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path peer 
    path gene_expr

    output:
    path "covariates.txt", emit: covariate_file
    
    script:
    """
    process_covariates.R ${peer} ${gene_expr} covariates.txt
    """
}

process generate_pcs {
    tag "Principal Components"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/Covariates" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input: 
    path gene_expr

    output:
    path "covariates.txt", emit: covariate_file
    
    script:
    """
    generate_pcs.R ${gene_expr} covariates.txt
    """
}
