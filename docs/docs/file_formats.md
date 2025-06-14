# File Formats

## Annotation Files
**remeta** uses many of the same files as **regenie** to define variant sets.
Most annotation files compatible with **regenie** should also be compatible with **remeta**.

### `--anno-file`
```
1:55039839:T:C PCSK9 LoF
1:55039842:G:A PCSK9 missense
.
```
A file defining variant annotations.
Contains 3 whitespace delimited columns: variant id (in CPRA format), gene name, and variant annotation.

**New in v0.9.0**: **remeta** now supports a tabixible 5-column annotation file where column 4 is the chromosome and column 5 is the position.
```
1:55039839:T:C PCSK9 LoF 1 55039839
1:55039842:G:A PCSK9 missense 1 55039842
.
```

### `--set-list`
```
A1BG 19  58346922  19:58346922:C:A,19:58346924:G:A,...
A1CF 10  50806630  10:50806630:A:G,10:50806630:A:AT,...
.
```
A file defining which variants are part of a gene set.
Contains 4 whitespace delimited columns: gene name, chromosome, start position, and a comma separated list of variants in the gene.

### `--mask-def`
```
LoF LoF
missense missense
LoF_missense LoF,missense
.
```
A file specifying which annotation categories to combine into masks.
Contains 2 whitespace delimited columns: mask name, and a comma separated list of variant annotations in a mask.

### `--aaf-file`
An optional file with variant alternate allele frequencies.
If specified, these frequencies are used for building masks.
Three formats are supported.
All formats assume whitespace separated columns.

* **Two Column Format**: CPRA and allele frequency
* **Five Column Format**: CPRA, frequency, is singleton (either 0 or 1), chromosome, position
* **Six Column Format** (plink2 afreq file): chromosome, CPRA, ref, alt, alt_freqs, and obs_ct

Here is an example of the five column format:
```
7:6187101:C:T 1.53918207864341e-05 0 7 6187101
7:6190395:C:A 2.19920388819247e-06 1 7 6190395
.
```

## Reference LD Files

### `--ld-prefix`
A set of three files named `$PREFIX.remeta.gene.ld`, `$PREFIX.remeta.buffer.ld`, and `$PREFIX.remeta.ld.idx.gz` generated by `remeta compute-ref-ld`.
The index `$PREFIX.remeta.ld.idx.gz` is bgzipped and human readable.
The columns are:

* `GENE_NAME`
* `GENE_LD_INDEX`
* `BUFFER_LD_INDEX`
* `GENE_VARIANTS_STORED`
* `BUFFER_VARIANTS_STORED`

### `--gene-list`
```
PCSK9	1	55039446	55064852
USP24	1	55066358	55215364
.
```
A file listing gene start and end positions.
Contains 4 whitespace separated columns: gene name, chromosome, start position, end position.
Note that this file must align with the `--set-list` file for gene-based tests.
Specifically,

1. Gene names in the `--gene-list` file should match gene names in the `--set-list` file exactly.
2. Any variant in the `--set-list` must appear within the start position and end position of the LD matrix for that gene: otherwise it will be dropped from the test (unless the `--keep-variants-not-in-ld-mat` option is specified).
3. The start and end position should be inclusive of any variant that could appear in a test.
In particular, it is not recommended to set the start and end positions based on the variants in the set list.
If the set list changes, this could result in variants being dropped from a test.

In the **remeta** repository we provide an example gene list under `resources/Ensembl100.GRCh38.chr1_23.gene_list.txt.gz` directory.
This file was created by extracting gene boundaries for protein coding genes annotated in Ensembl release 100.
An Ensembl GTF file [was downloaded from Ensembl](https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz);
start and end positions were extracted from lines with feature `gene` and the biotype `protein_coding`.

### `--target-pfile` and `--buffer-pfile`
A set of [pgen, pvar, and psam](https://www.cog-genomics.org/plink/2.0/input#pgen) files from plink2.

### `--genetic-map`
A genetic map in the SHAPEIT format.
Contains 3 columns: position, chromosome, and centimorgan.
Note that genetic maps are available from the [SHAPEIT5 repository](https://github.com/odelaneau/shapeit5).

## Miscellaneous Files

### `--htp`
A file in `htp` format, the default output format of **remeta**.
Output by **regenie** with the `--htp` option.
`htp` files are whitespace separated with the following columns:

* `Name`: Variant name. For **remeta**, single variants must be encoded `CHROM:POS:REF:ALT`, matching the naming convention of the annotation files.
* `Chr`: Chromosome.
* `Pos`: Position.
* `Ref`: Reference allele.
* `Alt`: Alternate allele.
* `Cohort`: Name of the cohort (e.g. UK Biobank).
* `Model`: Association model.
* `Effect`: Effect size. Note that for binary traits these are odds-ratios, not log odds ratios.
* `LCI_effect`: Lower 95% confidence interval for the effect size.
* `UCI_effect`: Upper 95% confidence interval for the effect size.
* `Pval`: Association p-value.
* `AAF`: Alternate allele frequency.
* `Num_Cases`: Number of cases. For QTs this is the sample size.
* `Cases_Ref`: Number of cases homozygous for the reference allele.
* `Cases_Het`: Number of cases who are heterozygotes.
* `Cases_Alt`: Number of cases homozygous for the alternate allele.
* `Num_Controls`: Number of controls. For QTs this is NA (as are the fields below).
* `Controls_Ref`: Number of controls homozygous for the reference allele.
* `Controls_Het`: Number of controls who are heterozygotes.
* `Controls_Alt`: Number of controls who are homozygous for the alternate allele.
* `Info`: A semi-colon (;) separated list of arbitrary fields. These can be `KEY=VALUE` pairs or a single field alone.

### `--extract` and `--exclude`
Files with variant ids (one per line) to include or exclude from meta-analysis.

### `--genep-def`
```
LoF LoF
NonSyn LoF,missense,LoF_missense
.
```
A file defining which masks to combine with ACAT.
Contains two space-separated columns: the name of the GENEP set, and a comma separated list of masks to include.

