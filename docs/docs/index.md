# remeta

**remeta** is a collection of C++ programs for meta-analysis of single variants and gene-based tests using summary statistics from [**regenie**](https://rgcgithub.github.io/regenie/).
It is primarily designed to perform meta-analysis of gene-based tests using exome variants,
and to perform conditional analysis of gene-based tests using imputed variants.

The main features of **remeta** are:

* It performs burden, ACATV, and SKATO meta-analysis from single variant summary statistics.
* Annotations and gene-sets can be updated without re-analysis of individual-level data.
* It only needs a single LD file per study.
* It provides an efficient utility to compute and store LD matrices.
* It can meta-analyze both quantitative and binary traits.
* It very accurately estimates effect sizes, allele frequencies, and genotype counts of gene burden masks from single variant summary statistics.
* It performs effect-size meta-analysis and p-value meta-analysis of single variants.

See the [**remeta** tutorial](tutorial.md) for a step-by-step example.

## Citation
Joseph, T., Mbatchou, J., et al. Computationally efficient meta-analysis of gene-based tests using summary statistics in large-scale genetic studies. Nat Genet (2025). [https://doi.org/10.1038/s41588-025-02390-0](https://doi.org/10.1038/s41588-025-02390-0).

## License
**remeta** is distributed under an MIT license.

## Contact
If you have any question about **remeta** please contact

* <jonathan.marchini@regeneron.com>
* <tyler.joseph@regeneron.com>

Please submit issues about **remeta** on [Github](https://github.com/rgcgithub/remeta/issues).