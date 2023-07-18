# PredictDb-nextflow: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
  * [`--gene_annotation`](#--gene_annotation)
  * [`--snp_annotation`](#--snp_annotation)
  * [`--genotype`](#--genotype)
  * [`--gene_exp`](#--gene_exp)
* [Other command line parameters](#other-command-line-parameters)
  * [`--covariates`](#--covariates)
  * [`--pca/--peer`](#--pca/--peer)
  * [`--outdir`](#--outdir)
  * [`--keepIntermediate`](#--keepIntermediate)
  * [`--prefix`](#--prefix)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO nf-core: Document required command line parameters to run the pipeline-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run main.nf --gene_annotation 'gene_annot.gtf' --snp_annot 'snp_annnotation_file.vcf' --genotype 'genotype_file' --gene_exp 'Normalized_gene expression.csv'
```

This will launch the pipeline with using the local executor.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Main arguments

### `--gene_annotation`

Use this to specify the location of your input gene annotation file. For example:

```bash
--reads 'path/to/data/gene_annotation.gtf'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The gene annotation file should contain annotations for all genes in 22 chromosomes
3. The gene annotation file should follow this format described [here]()

### `--snp_annotation`

Use this to specify the location of your input SNP annotation file. For example:

```bash
--snp_annot 'path/to/snp_annotation.vcf'
```

1. The path must be enclosed in quotes
2. The SNP annotation file should contain annotations for all genes in 22 chromosomes
3. The SNP annotation file should follow this format described [here](https://github.com/hakyimlab/PredictDb-nextflow/blob/master/docs/file_descriptions.md#snp-annotation-file)

### `--genotype`

Use this to specify the location of your genotype/dosage file. For example:

```bash
--genotype 'path/to/genotype_file.txt'
```
1. The genotype file contains the dosage of each sample for the specific varID and must be provided
2. It should have samples on the columns and the varID on the rows
3. Further description can be found [here](https://github.com/hakyimlab/PredictDb-nextflow/blob/master/docs/file_descriptions.md#genotype-file)

### `--gene_exp`

Use this to specify your gene expression file. For example:

```bash
--gene_exp 'path/to/gene_expression_file.txt'
```
1. The gene expression file must be normalized
2. Samples should  be on the columns while TargetID are on the rows
3. More indepth description of the gene expression file and preprocesing can be found [here](https://github.com/hakyimlab/PredictDb-nextflow/blob/master/docs/file_descriptions.md#gene-expression-file)


## Other command line parameters

### `--covariates`
The covariates to be regressed out from the gene expression

### `--pca/--peer`
Compute the principal components or peer factors to be regressed out from the the gene expressiom. You can either `--pca` or `--peer` not both.
This can be used in combination with `--covariates`. The `--pca` by default uses the first 10 principal components.

### `--outdir`

The output directory where the results will be saved.

### `--keepIntermediate`

By default the execution doesn't give you the intermediate files only the final output of the workflow. If you want to have intermediate files in your `outdir` provide this parameter.

### `--prefix`
Use this command to input the prefix name of your output files of the trained models. We recommend this to be the name of the population you are training the model on. For example:

```bash
--prefix 'Europeans'
```
If not provided it will use the provided generic name and will overwrite the initial output if the outdir is not provided.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

**NB:** Ensure mail or sendmail is set up correctly in your host before using this argument.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.


