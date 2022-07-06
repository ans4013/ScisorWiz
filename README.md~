ScisorWiz: Vizualizer for Differential Isoform Expression
================

## README

ScisorWiz is a linux-based R-package for visualizing differential
isoform expression for any gene across any number of cell types.

ScisorWiz affords its user the opportunity to visualize specific genes
across multiple cell types, and provides various sorting options for the
user to appraise the data in different ways. ScisorWiz provides a clear
picture of cell-type specific isoform expression through various
clustering methods and highlighting features such as separate colors for
alternative exons and Single Nucleotide Variants (SNVs). Additionally,
this package includes a feature in which the user can choose to generate
an interactive plot for exploratory purposes. This will generate a .jpg
file version and an interactive html file version of the plot.

------------------------------------------------------------------------

## Hardware / software requirements

ScisorWiz is a lightweight package with minimal dependencies, however,
for the simplest use case we do require an installation of:

-   R >= 3.5
-   python >= 3.7

Depending on the user’s needs, there are a few optional steps that
require additional software installation:

-   for SNV calling and representing insertions, deletions, and
    substitutions, we require
    [samtools](https://github.com/samtools/samtools)
-   for an exploratory and interactive view of the isoform plot which
    allows for zooming and/or cropping functionality, we require the
    [plotly](https://plotly.com/r/) library along with a few other
    packages. If the interactive feature is turned on, ScisorWiz will
    automatically download these packages as dependencies.

This package has been optimized to run with a single input file which
contains all the information necessary for visualizing isoform
expression across cell types (AllInfo.gz). This file is produced by the
[scisorseqr](https://github.com/noush-joglekar/scisorseqr/) package, and
we recommend the installation for a hassle-free process. However, the
input for ScisorWiz can also be obtained by other means. For more
information, please consult the vignette.

The package has only been tested on a CentOS x86_64 machine.

## Installation

The easiest way to install ScisorWiz is through
[Github](https://github.com) with:

``` r
devtools::install_github('ans4013/ScisorWiz',build_vignettes = TRUE)
```

## Workflow

<img src="man/figures/ScisorWiz_Workflow.png" width="100%" />

Some example data is included with this package to allow the user to
explore the input data and reproduce the use case. For example, running
the **ScisorWiz_AllInfo** function can be running by providing a path to
the input file (AllInfo.gz), the gene name (“*Snap25*”), the clustering
parameter (1 for intron chain), and a mismatch file (optional)
containing SNV information (Snap25.mismatches.txt.gz) using the
following command:

``` r
library(ScisorWiz)


ScisorWiz_AllInfo(gencodeAnno = "gencode.vM21.annotation.gtf.gz", 
AllInfoInput = allInfoFile,
cellTypeFile = cTypeFile, 
gene = "Snap25", cluster = 1, ci = .05,
mismatchCutoff = .05,
outputDir = "outputDir/",
mismatchFile = mismatchesFile,
interactive = "n")
```

A step-by-step outline of the various functions is available as a
vignette. To access it, run

``` r
browseVignettes("ScisorWiz")
```

## Support

We appreciate any and all inputs for improving ScisorWiz. Feel free to
send us an [email](mailto:astein2050@gmail.com).
