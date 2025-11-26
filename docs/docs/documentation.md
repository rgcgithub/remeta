# Documentation
This pages documents the various options for each subprogram in **remeta**.

```none
Usage: remeta [OPTIONS] COMMAND [ARGS]...

Options:
  -h [ --help ]         Print this message and exit.
  -v [ --version ]      Print version.

Commands:

* Meta-Analyses
  pvma         Perform p-value meta-analysis.
  esma         Perform effect size meta-analysis.
  gene         Perform burden, SKATO, and ACATV meta-analysis.
  merge        Merge results and compute meta-analysis GENEP.

* Utilities
  compute-ref-ld   Compute a reference LD matrix from plink2 pgen/pvar/psam files.
  index-anno       Create an index from a bgzipped REGENIE annotation file.
```


## remeta gene
Perform burden[^1], SKATO[^2], and ACATV[^3] meta-analysis.

**Marginal testing**

Gene-based meta-analysis requires LD matrices computed by `remeta compute-ref-ld`.
Each matrix is represented by three files: `$PREFIX.remeta.gene.ld`, `$PREFIX.remeta.buffer.ld` `$PREFIX.remeta.ld.idx.gz`.
The prefix of each file is passed to the `--ld-prefixes` argument.
These must be in the same order as the `--htp` and `--cohorts` arguments.

```none
remeta gene \
  --htp HTP1 HTP2 ... \
  --ld-prefixes LD_FILE1 LD_FILE2 ... \
  --cohorts COHORT1 COHORT2 ... \
  --anno-file ANNO_FILE \
  --set-list SET_LIST_FILE \
  --mask-def MASK_DEF_FILE \
  --trait-name MY_TRAIT \
  --trait-type TYPE \
  --out OUT_PREFIX
```

**Conditional analysis**

Conditional analysis can be performed for gene-based tests using the `--condition-htp` and `--condition-list` arguments.
The `--condition-htp` argument takes one HTP file per cohort and can be the same files passed to `--htp`.
The `--condition-list` argument takes a file with variant IDs to condition on (one variant ID per line).
```none
remeta gene \
  --htp HTP1 HTP2 ... \
  --ld-prefixes LD_FILE1 LD_FILE2 ... \
  --cohorts COHORT1 COHORT2 ... \
  --anno-file ANNO_FILE \
  --set-list SET_LIST_FILE \
  --mask-def MASK_DEF_FILE \
  --trait-name MY_TRAIT \
  --trait-type TYPE \
  --out OUT_PREFIX \
  --condition-list VARIANT_ID_FILE \
  --condition-htp HTP1 HTP2 ...
```

**Running without LD matrices (not recommended)**

`remeta gene` can be run without the required LD matrices by specifying the `--ignore-mask-ld` and `--keep-variants-not-in-ld-mat` flags.
Note that it is not possible to perform conditional analysis without LD matrices.

```none
remeta gene \
  --htp HTP1 HTP2 ... \
  --ignore-mask-ld \
  --keep-variants-not-in-ld-mat \
  --cohorts COHORT1 COHORT2 ... \
  --anno-file ANNO_FILE \
  --set-list SET_LIST_FILE \
  --mask-def MASK_DEF_FILE \
  --trait-name MY_TRAIT \
  --trait-type TYPE \
  --out OUT_PREFIX
```

**Specifying allele frequencies**

**remeta** provides several options to specify allele frequencies used to build masks.
By default, **remeta** computes on overall allele frequency per variant based on all cohorts where the variant was observed.
Alternatively **remeta** can use a maximum allele frequency observed across cohorts with the `--af-strategy max` argument.
Lastly, **remeta** allele frequencies can be specifed in an allele frequency file using the `--aaf-file` argument.
See [File Formats](file_formats.md) for a list of available formats.

**Variant weights**

Weights can be adjusted with the `--{skato,burden}-weight-strategy` argument. Setting `--{skato,burden}-weight-strategy beta` draws weights $w^{1/2} = \text{Beta}(MAF; 1, 25)$. Setting `--{skato,burden}-weight-strategy uniform` sets all variant weights to 1.

**Unbalanced binary traits**

**remeta** uses two strategies to control type 1 error for unbalanced binary traits:
a saddlepoint approximation (SPA) applied per mask and an SPA applied per variant.
The default parameters apply a mask level or variant level SPA when the case-control ratio of the trait falls below a certain threshold.
Simulations suggest that the threshold on case control apply an SPA depends on the test (e.g. burden vs. SKATO),
so parameters can be adjusted per test using several command line parameters.
Mask level parameters are available for burden tests and SKATO, and adjusted using the `--<burden,skato>-mask-spa-<pval,ccr>` arguments.
Variant level parameters are avaiable for burden test, SKATO, and ACATV, and adjusting using the `--<burden,skato,acatv>-sv-spa-<pval,ccr>` arguments.

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--htp`                     | FILE1 FILE2 ... | Required | List of HTP input files. |
| `--ld-prefixes`             | FILE1 FILE2 ... | Required | Prefix to `.remeta.gene.ld`, `.remeta.buffer.ld`, and `.remeta.ld.idx.gz` files per cohort. |
| `--cohorts`                 | STRING1 STRING2 ... | Required | List of cohort names per file. |
| `--anno-file`               | FILE | Required | File with variant annotations. Bgzipped and and indexed with `index-anno`. |
| `--set-list`                | FILE | Required | Regenie set-list file. |
| `--mask-def`                | FILE | Required | Regenie mask-def file. |
| `--trait-name`              | STRING | Required | Name of trait. |
| `--trait-type`              | STRING | Required | One of BT or QT. |
| `--out`                     | STRING | Required | Prefix for output files. |
| `--burden-aaf-bins (=0.0001 0.001 0.005 0.01)` | FLOAT1 FLOAT2 ... | Optional | Allele frequency cutoffs for building masks for burden testing. |
| `--burden-singleton-def (=within)` | STRING | Optional | Define singletons for the singleton mask within cohorts or across cohorts. One of 'within', 'across' or 'omit'. |
| `--burden-weight-strategy (=uniform)` | STRING | Optional | Strategy to compute variant weights for burden testing. One of `beta` or `uniform`. |
| `--burden-mask-spa-pval (=0.05)` | FLOAT | Optional | Apply a mask level SPA to burden tests when p-value < spa pval (BTs only). |
| `--burden-mask-spa-ccr (=0.01)`  | FLOAT | Optional |  Apply a mask level SPA to burden tests # cases / # controls < spa-ccr (BTs only). |
| `--burden-sv-spa-pval (=0.05)`   | FLOAT | Optional |  Apply a per variant SPA to burden tests p-value <  spa pval (BTs only). |
| `--burden-sv-spa-ccr (=0.00)`    | FLOAT | Optional |  Apply a per variant SPA to burden tests when # cases / # controls < spa ccr (BTs only). |
| `--skip-burden`              | FLAG | Optional | Do not run burden testing. |
| `--skato-max-aaf (=0.01)`    | FLOAT | Optional | Maximum allele frequency for a variant to be included in mask for SKATO. |
| `--skato-rho-values (=0 0.01 0.04 0.09 0.16 0.25 0.5 1)` | FLOAT1 FLOAT2 ... | Optional | Rho values for SKATO. |
| `--skato-min-aac (=1)`      | INT | Optional |  Minimum AAC across cohorts for a variant to be included in a mask for SKATO. |
| `--skato-weight-strategy`   | STRING | Optional | Strategy to compute variant weights for SKATO. One of 'beta' or 'uniform'. |
| `--skato-mask-spa-pval (=0.05)` | FLOAT | Optional |  Apply a mask level SPA to SKATO when p-value < spa pval (BTs only). |
| `--skato-mask-spa-ccr (=0.02)`  | FLOAT | Optional |  Apply a mask level SPA to SKATO when # cases / # controls < spa ccr (BTs only). |
| `--skato-sv-spa-pval (=0.05)`   | FLOAT | Optional |  Apply a per variant SPA to SKATO when p-value < spa pval (BTs only). |
| `--skato-sv-spa-ccr (=0.02)`    | FLOAT | Optional |  Apply a per variant SPA to SKATO when #cases / # controls < spa ccr (BTs only). |
| `--skip-skato`              | FLAG | Optional | Do not run SKATO. |
| `--acatv-max-aaf (=0.01)`   | FLOAT | Optional | Maximum allele frequency for a variant to be included in mask for ACATV. |
| `--acatv-min-aac (=5)`      | INT | Optional | Minimum AAC across cohorts for a variant to be included in a mask for ACATV. |
| `--acatv-weight-strategy`   | STRING | Optional | Strategy to compute variant weights for ACATV. One of 'beta' or 'uniform'. |
| `--acatv-sv-spa-pval (=0.05)` | FLOAT | Optional | Apply a per variant SPA to ACATV when p-value < spa pval (BTs only). |
| `--acatv-sv-spa-ccr (=0.02)`  | FLOAT | Optional | Apply a per variant SPA to ACATV when #cases / # controls < spa- ccr (BTs only). |
| `--skip-acatv`              | STRING | Optional | Do not run ACATV. |
| `--condition-list`          | FILE | Optional | File with variants to condition on (one per line). |
| `--condition-htp`           | FILE1 FILE2 ... | Optional | List of HTP files with summary statistics of conditional variants per cohort. |
| `--af-strategy (=overall)`  | STRING | Optional |  Strategy to compute variant allele frequences. One of 'overall' or 'max'. |
| `--aaf-file`                | FILE | Optional | Use precomputed alternate allele frequencies from an external file. |
| `--chr`                     | STRING | Optional | Run only on specifed chromosome. |
| `--gene`                    | STRING | Optional | Run only on specified gene. |
| `--extract`                 | FILE | Optional | Include only the variants with IDs listed in this file (one per line). |
| `--cohort-extract`          | COHORT1=FILE1 COHORT2=FILE2 ... | Optional | Only include variants in FILE_N from COHORT_N and not in --exclude. |
| `--exclude`                 | FILE | Optional | Exclude variants with IDs listed in this file (one per line, overridden by --extract). |
| `--sources`                 | STRING1 STRING2 ... | Optional | Only meta-analyze variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...|
| `--write-cohort-burden-tests` | FLAG | Optional |  Compute and store per cohort burden tests (ignores changes to --burden-weight-strategy). |
| `--write-mask-snplist`      | FLAG | Optional | Write file with list of variants included in each mask. |
| `--recompute-score`         | FLAG | Optional | Recompute score statistics from betas and standard errors when missing in input. |
| `--keep-variants-not-in-ld-mat` | FLAG | Optional | Keep variants absent from the LD matrix instead of dropping them. | 
| `--ignore-mask-ld`          | FLAG | Optional | Ignore LD between variants in a mask. |
| `--write-variant-aaf`       | FLAG | Optional | Output variant AAFs used to construct masks. |
| `--threads (=1)`            | INT | Optional | Number of threads to use. |



## remeta esma
Perform effect size meta-analysis (i.e. inverse variance meta-analysis).

**Effect-size meta-analysis of single variants**
```none
remeta esma \
  --htp HTP1 HTP2 ... \
  --cohorts COHORT1 COHORT1 ... \
  --trait-name MY_TRAIT \
  --trait-type TYPE \
  --out OUT_PREFIX
```

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--htp`                     | FILE1 FILE2 ... | Required | List of HTP input files. |
| `--cohorts`                 | STRING1 STRING2 ... | Required | List of cohort names per file. |
| `--trait-name`              | STRING | Required | Name of trait. |
| `--trait-type`              | STRING | Required | One of BT or QT. |
| `--out`                     | STRING | Required | Prefix for output files. |
| `--chr`                     | STRING | Optional | Run only on chromosome CHR. |
| `--extract`                 | FILE | Optional | Include only the variants with IDs listed in this file (one per line). |
| `--exclude`                 | FILE | Optional | Exclude variants with IDs listed in this file (one per line). |
| `--sources`                 | STRING1 STRING2 ... | Optional | Only meta-analyze variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...|
| `--source-def`              | FILE | Optional | Two column file mapping long SOURCE info fields to shorter SOURCE info fields. |


## remeta pvma
Perform p-value meta-analysis with either Stouffer's [^4] method or Fisher's method [^5].
Does not currently use effect direction.
Primary use case is to meta-analyze ACATV, SKATO, and SBAT from **regenie** using standard meta-analysis.

It is recommended to run this command with the `--skip-beta` flag to avoid meta-analyzing single variants with both `pvma` and  `esma`.

**P-value meta-analysis of ACATV, SKATO, and SBAT (if ran in regenie)**
```none
remeta pvma \
  --htp HTP1 HTP2 ... \
  --cohorts COHORT1 COHORT1 ... \
  --trait-name MY_TRAIT \
  --trait-type TYPE \
  --out OUT_PREFIX \
  --skip-beta \
  --method stouffers
```

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--htp`                     | FILE1 FILE2 ... | Required | List of HTP input files. |
| `--cohorts`                 | STRING1 STRING2 ... | Required | List of cohort names per file. |
| `--trait-name`              | STRING | Required | Name of trait. |
| `--trait-type`              | STRING | Required | One of BT or QT. |
| `--out`                     | STRING | Required | Prefix for output files. |
| `--method (=stouffers)`     | STRING | Optional | One of `stouffers` or `fishers`. |
| `--unweighted`              | FLAG | Optional | Omit sample size weighting (affects stouffers only). |
| `--two-sided`               | FLAG | Optional | Use two-sided p-values when variants have an effect size in each cohort (affects stouffers only). |
| `--chr`                     | STRING | Optional | Run only on chromosome CHR. |
| `--extract`                 | FILE | Optional | Include only the variants with IDs listed in this file (one per line). |
| `--exclude`                 | FILE | Optional | Exclude variants with IDs listed in this file (one per line). |
| `--skip-beta`               | FLAG | Optional | Skip entries with an effect size estimate. |
| `--sources`                 | STRING1 STRING2 ... | Optional | Only meta-analyze variants where the info field SOURCE is one of SOURCE1 SOURCE2 ...|
| `--source-def`              | FILE | Optional | Two column file mapping long SOURCE info fields to shorter SOURCE info fields. |


## remeta merge
Merge results and compute meta-analysis GENE_P [^6].

**Compute GENE_P from remeta's gene-based tests**
```none
remeta merge \
  --htp PVMA_HTP ESMA_HTP GENE_HTP1 ... GENE_HTP23  \
  --genep-def GENEP_DEF_FILE \
  --out OUT_PREFIX
```

**Compute GENE_P from regenie's gene-based tests (QT with additive model)**
```none
remeta merge \
  --htp GENE_HTP1 ... GENE_HTP23  \
  --genep-def GENEP_DEF_FILE \
  --burden-model ADD-WGR-LR \
  --acatv-model ADD-WGR-ACATV \
  --skato-model ADD-WGR-SKATO-ACAT \
  --sbat-model ADD-WGR-BURDEN-SBAT \
  --out OUT_PREFIX
```

**Compute GENE_P from regenie's gene-based test (BT with additive model using Firth regression)**
```none
remeta merge \
  --htp GENE_HTP1 ... GENE_HTP23  \
  --genep-def GENEP_DEF_FILE \
  --burden-model ADD-WGR-FIRTH \
  --acatv-model ADD-WGR-ACATV \
  --skato-model ADD-WGR-SKATO-ACAT \
  --sbat-model ADD-WGR-BURDEN-SBAT \
  --out OUT_PREFIX
```

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--htp`                     | FILE1 FILE2 ... | Required | List of HTP input files. |
| `--out`                     | STRING | Required | Prefix for output files. |
| `--genep-def`               | STRING | Optional | File with masks to group for GENE_P. |
| `--chr`                     | STRING | Optional | Run only on specified chromosome. |
| `--burden-model (=REMETA-BURDEN-META)` | STRING | Model column to collapse burden p-values. |
| `--acatv-model (=REMETA-ACATV-META)` | STRING | Model column to collapse ACATV p-values. |
| `--skato-model (=REMETA-SKATO-META)` | STRING | Model column to collapse SKATO p-values. |
| `--sbat-model (=ADD-WGR-BURDEN-SBAT)` | STRING | Optional | Include SBAT PVMA in GENE_P if available. |


## remeta compute-ref-ld
Compute reference LD matrices from plink2 pgen/pvar/psam files.

**Reference LD for marginal testing**
```none
remeta compute-ref-ld \
  --target-pfile PFILE \
  --gene-list GENE_LIST_FILE \
  --skip-buffer \
  --out OUT_PREFIX
```

**Reference LD matrices for conditional analysis**

If imputed variants are in a separate file from the WES variants:
```none
remeta compute-ref-ld \
  --target-pfile WES_PFILE \
  --buffer-pfile IMPUTED_PFILE \
  --gene-list GENE_LIST_FILE \
  --buffer-mb 1 \
  --out OUT_PREFIX
```

If imputed variants are in the same file as the WES variants:
```none
remeta compute-ref-ld \
  --target-pfile COMBINED_PFILE \
  --target-extract WES_VARIANT_LIST \
  --gene-list GENE_LIST_FILE \
  --buffer-mb 1 \
  --out OUT_PREFIX
```

Note that buffer regions can also be defined in centimorgans using the `--buffer-cm` option. This requires the `--genetic-map` argument with a genetic map file in the [SHAPEIT format](https://github.com/odelaneau/shapeit5/tree/main/resources/maps/b38).

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--target-pfile`            | STRING | Required | Prefix to pgen/pvar/psam files in target regions (typically WES).
| `--gene-list`               | FILE | Required | List of genes to include in LD matrix in four column format: GENE_NAME CHR START END. |
| `--chr`                     | STRING | Required | Chromosome to run. |
| `--out`                     | STRING | Required | Prefix for output files. |
| `--buffer-pfile`            | STRING | Optional | Prefix to pgen/pvar/psam files to use for buffer regions (typically array or imputed genotypes). |
| `--buffer-mb`               | FLOAT | Optional | Buffer in Mb around each gene to search for variants to include LD file. |
| `--buffer-cm`               | FLOAT | Optional | Buffer in cM around each gene to search for variants to include in LD file (requires a genetic map). |
| `--genetic-map`             | FILE | Optional | Path to genetic map in three column format: POS CHR CM. Required for `--buffer-cm` option. |
| `--target-r2 (=0.0001)`     | FLOAT | Optional | Drop target (gene) LD matrix entries where r2 < target_r2. |
| `--buffer-r2 (=0.0001)`     | FLOAT | Optional | Drop buffer (conditional) LD matrix entries where r2 < buffer_r2. |
| `--float-size (=1)`         | INT | Optional | Size of float in bytes to store LD of buffer variants. Possible values: 1, 2, or 4 |
| `--block-size (=2048)`      | INT | Optional | Number of genotypes loaded into memory = 2*block_size. |
| `--threads (=1)`            | INT | Optional | Number of threads for computation. | 
| `--target-extract`          | FILE | Optional | Extract list of variants to include in target regions (e.g. exonic variants) |
| `--target-exclude`          | FILE | Optional | Exclude list of variants from target file (e.g. non-coding variants) |
| `--target-keep`             | FILE | Optional | File with list of samples to keep (one sample per line, matching columns of psam file) |
| `--target-remove`           | FILE | Optional | File with list of samples to remove (one per line, matching columns of psam file) |
| `--buffer-extract`          | FILE | Optional | Extract list of variants to include in buffer regions (e.g. imputed variants) | 
| `--buffer-exclude`          | FILE | Optional | Exclude list of variants from buffer file (e.g. low quality variants) |
| `--skip-buffer`             | FLAG | Optional | Exclude all buffer variants from LD calculation. |


## remeta index-anno
Index annotation files. Files must be bgzipped

```none
remeta index-anno --file ANNOTATION_FILE.gz
```

### Options

| Option | Argument | Type | Description |
|--------|----------|------|-------------|
| `--file`                    | FILE | Required | Path to annotation file. |
| `--stride`                  | INT | Required | Length in bases between each index pointer. |

## References

[^1]: Lee, S., Abecasis, G. R., Boehnke, M., & Lin, X. (2014). Rare-variant association analysis: study designs and statistical tests. The American Journal of Human Genetics, 95(1), 5-23.

[^2]: Lee, S., Emond, M. J., Bamshad, M. J., Barnes, K. C., Rieder, M. J., Nickerson, D. A., ... & Lin, X. (2012). Optimal unified approach for rare-variant association testing with application to small-sample case-control whole-exome sequencing studies. The American Journal of Human Genetics, 91(2), 224-237.

[^3]: Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin, X. (2019). ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. The American Journal of Human Genetics, 104(3), 410-421.

[^4]: [Stouffer's Method](https://en.wikipedia.org/wiki/Fisher%27s_method#Relation_to_Stouffer's_Z-score_method)

[^5]: [Fisher's Method](https://en.wikipedia.org/wiki/Fisher%27s_method)

[^6]: Ziyatdinov, A., Mbatchou, J., Marcketta, A., Backman, J., Gaynor, S., Zou, Y., ... & Marchini, J. (2024). Joint testing of rare variant burden scores using non-negative least squares. The American Journal of Human Genetics, 111(10), 2139-2149.