## Input files
1. **Gene Expression file** - A dataframe/matrix where typically rows are genes and columns are subjects with values being inverse quantile normalized expression values. 
This fileusually comes in tab delimited text file or .RDS files.
2. **Genotype file** - A matrix where rows are variant ID’s of SNPs and columns are subjects with values between 0 and 2 indicating the dosage for the second allele listed in the 
variant ID. This could come in a .vcf, .RDS, or .txt file format.
3. **SNP annotation file** - A dataframe where rows are SNPs, columns are: chromosome,position, VariantID (chr_pos_refAllele_altAllele_build) RefAllele, AltAllele, rsid
4. **Gene annotation file** - A dataframe containing information about each gene. Rows are genes with ensemble ID’s, columns are chr, gene_id, gene_name, start, end. 
This file is most likely a .gtf file but could also be a .RDS file.
5. **Covariate file** - A matrix with the covariates to regress out of the gene expression matrix. Each row is a subject and the columns are covariates such as sex, age, 
and PEER factors. This is usually a tab delimited text file but most of the times we will need to create this file ourselves. Its purpose is to run regressions for each gene in the gene expression file and extract the residuals of that regression to regress out the effects of the covariates for each 
(gene, sample) pair.

## General Preprocessing
Essentially, we want to get all of our input files from whatever format they come to us into the way that our pipeline takes it. 
We want each file to look like this:
- **Gene Annotation:** an R dataframe saved as an RDS object with columns chr, gene_id, gene_name,start, end. The row names of the dataframe are the gene_id.
- **SNP Annotation:** The snp annotation file is split up by chromosome (i.e. there are 22 files. There's an assumption of working with human genotypes,
and we don't work with X or Y chromosome and saved as an R dataframe (RDS Object) with columns chr, pos,varID, refAllele, effectAllele, and rsid. 
All rows are SNPs (no INDELS), andunambiguously stranded. The varID column is saved as the rownames.
- **Genotype files:** again split by chromosome, but saved as a tab-delimited text file. Has a header row with fields Id for the variant ID of the snp and the rest are the 
individual ID labels. Rows have been filtered so that only single letter variants are included.
- **Expression files:** An R dataframe saved as an RDS object. colnames are the ensembl ID and the row names are the individual ids. The values are normalized with covariate 
effects regressed out. *__Note__* that this is transposed from the original input.
