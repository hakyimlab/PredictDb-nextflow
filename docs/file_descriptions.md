# Input files
### **Gene Expression file** 
A dataframe/matrix where typically rows are genes and columns are subjects with values being inverse quantile normalized expression values. 

Expression and covariates are tab-separated text files with gzip-compression. They look like:
```
NAME              IND1  IND2  IND3 ...

ENSG00000227232.5 -0.81 -0.29 0.31 ...
... 
```
The first row is the header. For expression, the first column must be NAME. Each row contains the gene expression values for the individuals in your sample. For covariates, the first column is ID and each row contains a covariate such as ancestral pca, peer factors, etc. The pipeline doesn't allow missing values and individuals must be the same in expression and covariate files.

### **Genotype file**
A matrix where rows are variant ID’s of SNPs and columns are subjects with values between 0 and 2 indicating the dosage for the second allele listed in the 
variant ID. This could come in a .vcf or .txt file format.

Measured genotypes are text files with the following format:
```
varID           IND1 IND2  IND3 ...

1_54421_A_G_b37 1    0     0    ... 

...
```
These files are similar to the expression and covariate files. Dosages are expected to be in the [0,2] range. Missing values are not supported.

They must be accompanied by variant annotation files (SNP annotation file). 

## **SNP annotation file**
A dataframe where rows are SNPs, columns are: chromosome,position, VariantID (chr_pos_refAllele_altAllele_build) RefAllele, AltAllele, rsid
They are text files with the following format:
```
chromosome pos    varID            ref_vcf alt_vcf R2                 MAF     rsid      rsid_dbSNP150

1          566875 1_566875_C_T_b37 C       T       0.9747600000000001 0.03085 rs2185539 rs2185539
...
```
The variant ids need not be in the "{chr}_{pos}_{NEA}_{EA}_b37" format; the variant annotation file must provide the mapping from variant id to rsid.

## **Gene annotation file** 
A dataframe containing information about each gene. Rows are genes with ensemble ID’s, columns are chr, gene_id, gene_name, start, end. 

This file is most likely a .gtf file or a text file containing the list of genes you want to train on. It can contain less entries than those available to the expression data files. The file looks like this:
```
chr gene_id           gene_name     start end    gene_type

1   ENSG00000243485.5 MIR1302-2HG   29554 31109  lincRNA

1   ENSG00000237613.2 FAM138A       34554 36081  lincRNA

1   ENSG00000186092.4 OR4F5         69091 70008  protein_coding

1   ENSG00000238009.6 RP11-34P13.7  89295 133723 lincRNA

...
```

## **Covariate file** 
A matrix with the covariates to regress out of the gene expression matrix. Each row is a subject and the columns are covariates such as sex, age, 
and PEER factors. 

This is usually a tab delimited text file but most of the times we will need to create this file ourselves. Its purpose is to run regressions for each gene in the gene expression file and extract the residuals of that regression to regress out the effects of the covariates for each 
(gene, sample) pair. 

The covariate file should have a similar structure as the expression file.

# General Preprocessing
Essentially, we want to get all of our input files from whatever format they come to us into the way that our pipeline takes it. 
We want each file to look like this:

The genotype and annotation files are expected to be split by chromosome.

- **Gene Annotation:** an R dataframe saved as an RDS object with columns chr, gene_id, gene_name,start, end. The row names of the dataframe are the gene_id.
- **SNP Annotation:** The snp annotation file is split up by chromosome (i.e. there are 22 files. There's an assumption of working with human genotypes,
and we don't work with X or Y chromosome and saved as an R dataframe (RDS Object) with columns chr, pos,varID, refAllele, effectAllele, and rsid. 
All rows are SNPs (no INDELS), andunambiguously stranded. The varID column is saved as the rownames.
- **Genotype files:** again split by chromosome, but saved as a tab-delimited text file. Has a header row with fields Id for the variant ID of the snp and the rest are the 
individual ID labels. Rows have been filtered so that only single letter variants are included.
- **Expression files:** An R dataframe saved as an RDS object. colnames are the ensembl ID and the row names are the individual ids. The values are normalized with covariate 
effects regressed out. *__Note__* that this is transposed from the original input.
