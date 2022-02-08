ScisorWiz: Vizualizer for Differential Isoform Expression
================

## README

ScisorWiz is a linux-based R-package for visualizing differential
isoform expression for any gene across up to six cell types.

ScisorWiz affords its user the opportunity to visualize specific genes
across specific cell types, and provides various sorting options for the
user to gain different ways to understand their data. ScisorWiz provides
a clear picture of differential isoform expression through various
clustering methods and highlighting features such as alternative exons
and Single Nucleotide Variants (SNVs).

------------------------------------------------------------------------

## Hardware / software requirements

The package has only been tested on a CentOS x86_64 machine. For
obtaining SNVs, deletions and insertions we require
[samtools](https://github.com/samtools/samtools) software installation.

To run the ScisorWiz_AllInfo function (recommended for faster and more
in-depth information and required for the different clustering methods)
we require [scisorseqr](https://github.com/noush-joglekar/scisorseqr/)
software installation to obtain the AllInfo file.

## Installation

The easiest way to install ScisorWiz is through
[Github](https://github.com) with:

``` r
devtools::install_github('ans4013/ScisorWiz',build_vignettes = TRUE)
```

## Workflow

<!-- ``` {r prettyGraph, echo=FALSE, out.width = '60%'} -->
<!-- knitr::include_graphics("man/figures/ScisorWiz_workflow.png") -->
<!-- ``` -->

These steps are available as functions in the package. For example,
ScisorWiz_AllInfo with clustering by intron chain (1) and mismatches can
be done using the following command:

``` r
library(ScisorWiz)
ScisorWiz_AllInfo(gencodeAnno = "gencode.vM21.annotation.gtf.gz", AllInfoInput = "AllInfo.gz",
cellTypeFile = "userInput/celltypeFile", gene = "Snap25", cluster = 1, ci = .05,
outputDir = "extdata/outputDir/", mismatchFile="extdata/outputDir/Snap25.mismatches.txt.gz")
```

A step-by-step outline of the various functions is available as a
vignette. To access it, run

``` r
browseVignettes("ScisorWiz")
```

## Support

We appreciate any and all inputs for improving ScisorWiz Feel free to
send us an [email](mailto:ans4013@med.cornell.edu).
