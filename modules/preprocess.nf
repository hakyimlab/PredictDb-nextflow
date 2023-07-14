#! /usr/bin/env nextflow


process gunzipFile {
    tag "$gz"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/unzipped-files" : params.outdir },
               saveAs: { params.keepIntermediate ? it : null }, mode: 'copy'
    input:
    path gz

    output:
    path "${gz.baseName}"

    script:
    """
    gunzip --verbose --stdout --force ${gz} > ${gz.baseName}
    """
}


process parse_geneAnnot {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/gtf" : false },    
	       saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path gtf

    output:
    path "gene_annot.parsed.txt", emit: parsed_annot

    script:
    """
    parse_gtf.py ${gtf} gene_annot.parsed.txt
    """
}

process split_snpAnnot {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/snp_annot" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path snp

    output:
    path "snp_annot.*.txt", emit: snp_files

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
    path genot

    output:
    path "genotype.*.txt", emit: genotype_files

    script:
    """
    split_genotype_by_chr.py ${genot} genotype
    """
}

process transpose_geneExpr {
    tag "pre-processing"
    publishDir path: { params.keepIntermediate ? "${params.outdir}/pre-proccessed/transpose_gene-expr" : false },
               saveAs: { params.keepIntermediate ? it : false }, mode: 'copy'
    input:
    path geneExp

    output:
    path("transposed_gene_exp.txt"), emit: tr_expr 
    path("transposed_gene_exp.csv"), emit: tr_expr_peer
    tuple stdout, emit: tr_expr_count

    script:
    """
    transpose_gene_expr.R ${geneExp} transposed_gene_exp.csv transposed_gene_exp.txt
    wc -l transposed_gene_exp.csv | awk '{\$1=\$1};1' | cut -d " " -f1 | tr -d \'\\n\'
    """
}