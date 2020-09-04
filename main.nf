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

    nextflow run main.nf --gtf 'gene_annot.vcf' --snp 'snp_annnotation_file' --genotype 'genotype_file' --gene_expr 'Normalized_gene expression'
    

    Mandatory arguments:
        --gtf [file]
        --snp [file]
        --genotype [file]
        --gene_exp [file]

    Options:
      --keepIntermediate [bool]              
      --population [str]
      --outdir [str]

                 

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
summary['Gene annotation']       = params.gtf
summary['SNP annotation']        = params.snp
summary['Expression file']       = params.gene_exp
summary['Genotype file']         = params.genotype
summary['Max Resources']         = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']            = params.outdir
summary['Launch dir']            = workflow.launchDir
summary['Working dir']           = workflow.workDir
summary['Script dir']            = workflow.projectDir
summary['User']                  = workflow.userName

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


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
    tuple stdout, path("transposed_gene_exp.csv") into tr_expr, gene_cols 

    script:
    """
    transpose_gene_expr.R ${geneExp} transposed_gene_exp.csv 
    wc -l transposed_gene_exp.csv | cut -d " " -f1 | tr -d \'\\n\'
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
    tuple val(tcount), file(csv:'transposed_file.csv') from tr_expr

    output:
    path "calculated_peers/X.csv" into peers

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

process linear_regression {
    tag "regression"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/Covariates" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path peer from peers
    tuple val(tcount), file(gene_expr:'transposed_file.csv') from gene_cols

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
    path "gtex_v7_${pop}.db" into filtering

    script:
    pop = params.prefix
    """
    make_db.R ${models} ${weights} ${chroms} ${pop}
    """
}


process filter_database {
    tag "database"
    publishDir path: "${params.outdir}/filtered_db",
               saveAs: it, mode: 'copy'

    input:
    path all_db from filtering

    output:
    path "gtex_v7_${pop}_filtered_signif.db" into filtered

    script:
    pop = params.prefix
    """
    filter_db.R ${all_db} gtex_v7_${pop}_filtered_signif.db
    """
}

// Get the chromosome number
def getChromID( file ) {
    file.name.toString().find(/(.chr)(\d+)(.txt)/) { match, pref, chrom, ext -> chrom }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[PredictDB-Nextflow] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[PredictDB-Nextflow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[PredictDb] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[PredictDb] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[PredictDb]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[PredictDb]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

