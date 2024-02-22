#! /usr/bin/env nextflow

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Get the chromosome number
def getChromID( file ) {
    file.name.toString().find(/(.chr)(\d+)(.txt)/) { match, pref, chrom, ext -> chrom }
}

//Import modules
include {gunzipFile as checkGTF; 
            gunzipFile as checkGENO;
            gunzipFile as checkExpr;
            gunzipFile as checkSNPAnnot;
            parse_geneAnnot; 
            split_snpAnnot;
            split_genotype; 
            transpose_geneExpr } from '../modules/preprocess.nf' addParams(outdir: "${params.outdir}")

include {model_training_w_covs; 
            model_training_wo_covs;
            model_training_w_covs_nested;
            model_training_wo_covs_nested } from '../modules/train_model.nf' addParams(outdir: "${params.outdir}")

include {collectModel_summaries;
            collectChrom_summaries;
            collectWeight_summaries;
            collectModel_covariances;
            make_database; 
            filter_database;
            filter_cov} from '../modules/create_db.nf' addParams(outdir: "${params.outdir}") 

include {generate_peer_factors;
            process_peer_factors;
            generate_pcs;
            combine_covs} from '../modules/covariates.nf' addParams(outdir: "${params.outdir}") 


workflow UNZIP {
   take:
	gtf
	snp_annot
	genotype
	geneExp

   main:
   if (hasExtension(params.gene_annotation, 'gz')) {
	gtf = checkGTF(gtf)}
   if (hasExtension(params.snp_annotation, 'gz')) {
	snp_annot = checkSNPAnnot(snp_annot)}
   if (hasExtension(params.gene_exp, 'gz')) {
	geneExp = checkExpr(geneExp)}
   if (hasExtension(params.genotype, 'gz')) {
	genotype = checkGENO(genotype)}

   emit:
        gtf = gtf
        snp_annot = snp_annot
        genotype = genotype
        geneExp = geneExp
}

workflow PREPROCESS {
    take:
        gtf
        snp_annot
        genotype
        geneExp
    main:
    parse_geneAnnot(gtf)
    split_snpAnnot(snp_annot)
    split_genotype(genotype)
    transpose_geneExpr(geneExp)

    emit:
        gene_annot = parse_geneAnnot.out.parsed_annot
        snp_files = split_snpAnnot.out.snp_files
        genotype_files = split_genotype.out.genotype_files
        tr_expr = transpose_geneExpr.out.tr_expr
        tr_expr_peer = transpose_geneExpr.out.tr_expr_peer
        tr_expr_count = transpose_geneExpr.out.tr_expr_count
}

workflow COVS {
    take:
        geneExp
        samp_size
    main:
	    if(params.peer) {
	        log.info "Using PEER to calculate covariates"
            generate_peer_factors(geneExp,samp_size)
            covs = process_peer_factors(generate_peer_factors.out.peers,geneExp)

        }
        else if (params.pca){
            log.info "Using PCA to calculate covariates"
            covs = generate_pcs(geneExp)
	    } else {
            log.info "No covariates calculate you can use --peer or --pca to capture underlying structure in expression"
            covs = null       
        }
    emit:
        covs = covs   
}

workflow COMBINE_COVS {
    take:
        computed_covs
        p_covs
    main:
        if(params.peer && params.covariates || params.pca && params.covariates) {
            println "Using both computed and provided covariates"
            covs = combine_covs(computed_covs,p_covs)
        }
        else if (params.covariates){
            println "Using the provided covariates"
            covs = p_covs
        } else if (params.peer || params.pca){
            println "Using the computed (peer/pca) covariates"
            covs = computed_covs
        } else {
            covs = null
            println "No covariates calculated or provided"
        }
    emit:
        covs = covs
}

workflow TRAIN_MODEL {
    take:    
        gene_annot
        snp_files
        genotype_files
        geneExp
        covs
    main:
        //split thre count and 

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

        //snp_genotype_files.view()

        if(params.peer || params.pca || params.covariates) {
            println "Run CV-Enet with covariates"
            if(params.nested_cv){
                println "Running nested CV-Enet"
                model = model_training_w_covs_nested(covs.first(),geneExp.first(),gene_annot.first(),snp_genotype_files)
            }else {
                model = model_training_w_covs(covs.first(),geneExp.first(),gene_annot.first(),snp_genotype_files)
            }
        } else {
            println "Run CV-Enet without covariates"
            if(params.nested_cv){
                println "Running nested CV-Enet"
                model = model_training_wo_covs_nested(geneExp.first(),gene_annot.first(),snp_genotype_files)
            } else {
                model = model_training_wo_covs(geneExp.first(),gene_annot.first(),snp_genotype_files)
            }
        }
    emit:
        chrom_summaries = model.chrom_summaries
        weight_summaries = model.weight_summaries
        model_summaries = model.model_summaries
        model_covariances = model.all_covariances  
}

workflow CREATE_DB {
    take:
        chrom_summaries
        weight_summaries
        model_summaries
        all_covariances

    main:
        collectModel_summaries(model_summaries)
        collectChrom_summaries(chrom_summaries)
        collectWeight_summaries(weight_summaries)
        make_database(collectModel_summaries.out.all_model_sum,
                      collectWeight_summaries.out.all_weight_sum,
                      collectChrom_summaries.out.all_chrom_sum)
        filter_database(make_database.out.whole_db)
        collectModel_covariances(all_covariances)
        filter_cov(filter_database.out.filtered_db,
                    collectModel_covariances.out.all_model_sum)
}
