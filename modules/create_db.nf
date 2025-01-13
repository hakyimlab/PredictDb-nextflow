#! /usr/bin/env nextflow

process collectModel_covariances {
    tag "database"
    label "process_high"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

    input:
    path covariance 

    output:
    path "predict_db_${pop}.txt.gz", emit: all_model_sum

    script:
    pop = params.prefix
    """
    covariance_summary.R $covariance
    mv Covariances.txt predict_db_${pop}.txt
    gzip predict_db_${pop}.txt
    """
}

process collectModel_summaries {
    tag "database"
    label "process_low"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

    input:
    path model 

    output:
    path "Model_summary.txt", emit: all_model_sum

    script:
    """
    model_summary.R $model
    """
}


process collectWeight_summaries {
    tag "database"
    label "process_medium"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

    input:
    path weight 

    output:
    path "Weight_summary.txt", emit: all_weight_sum

    script:
    """
    weight_summary.R $weight
    """
}

process collectChrom_summaries {
    tag "database"
    label "process_low"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

    input:
    path chrom 

    output:
    path "Chromosome_summary.txt", emit: all_chrom_sum

    script:
    """
    chrom_summary.R $chrom
    """
}

process make_database {
    tag "database"
    label "process_medium"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/database" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'

    input:
    path models 
    path weights 
    path chroms 

    output:
    path "predict_db_${pop}.db", emit: whole_db

    script:
    pop = params.prefix
    """
    make_db.R ${models} ${weights} ${chroms} ${pop}
    """
}


process filter_database {
    tag "database"
    label "process_medium"
    publishDir path: "${params.outdir}/filtered_db",
               saveAs: it, mode: 'copy'

    input:
    path all_db

    output:
    path "predict_db_${pop}_filtered.db", emit: filtered_db

    script:
    pop = params.prefix
    """
    filter_db.R ${all_db} predict_db_${pop}_filtered.db
    """
}


process filter_cov {
    tag "database"
    label "process_medium"
    publishDir path: "${params.outdir}/filtered_db",
               saveAs: it, mode: 'copy'

    input:
    path filtered_db
    path all_covs

    output:
    path "predict_db_${pop}_filtered.txt.gz", emit: filtered_covs

    script:
    pop = params.prefix
    """
    filter_cov.R ${filtered_db} ${all_covs} predict_db_${pop}_filtered.txt
    gzip predict_db_${pop}_filtered.txt
    """
}
