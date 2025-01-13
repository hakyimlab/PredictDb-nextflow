#! /usr/bin/env nextflow

process generate_peer_factors {
    tag "PEER Factors"
    label "process_medium"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/PEER" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    input: 
    path csv
    val tcount

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
    label "process_low"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/covariates" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    input:
    path peer 
    path gene_expr

    output:
    path "peer_covariates.txt", emit: covariate_file
    
    script:
    """
    process_covariates.R ${peer} ${gene_expr} peer_covariates.txt
    """
}

process generate_pcs {
    tag "Principal Components"
    label "process_medium"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/covariates" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    input: 
    path gene_expr

    output:
    path "pca_covariates.txt", emit: covariate_file
    
    script:
    """
    generate_pcs.R ${gene_expr} pca_covariates.txt
    """
}

process combine_covs {
    tag "Combining all covariates"
    label "process_low"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/covariates" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    input: 
    path computed_covs
    path provided_covs

    output:
    path "all_covariates.txt", emit: covariate_file
    
    script:
    """
    combine_covariates.R ${computed_covs} ${provided_covs}  all_covariates.txt
    """
}
