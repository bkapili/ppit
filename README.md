
<!-- README.md is generated from README.Rmd. Please edit that file -->

![PPIT
logo](https://raw.githubusercontent.com/BKapili/ppit/master/logo/PPIT_logo_v1.png)

<!-- badges: start -->

<!-- badges: end -->

Phylogenetic Placement for Inferring Taxonomy (PPIT) is an R package for
inferring microbial taxonomy from *nifH* amplicons. PPIT provides a
higher proportion of correct taxonomic inferences at all taxonomic ranks
in comparison to BLAST-based approaches at the cost of fewer total
inferences.

## Installation

To install `ppit` directly from RStudio, first download the dependent
package `Biostrings` from Bioconductor then use the `devtools` package
to install `ppit`:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("BKapili/ppit")
```

## Documentation

To get started, see the
[tutorial](https://github.com/BKapili/ppit/blob/master/tutorial/ppit_tutorial.md)
that walks through the analysis using *nifH* amplicons generated from
Illumina MiSeq data.
