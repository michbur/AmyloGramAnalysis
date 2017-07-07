# Read me

This repository contains data and code necessary to reproduce analysis from the article: Burdukiewicz M, Sobczyk P, Rödiger S, Duda-Madej A, Mackiewicz P, Kotulska M. (2017) Amyloidogenic motifs revealed by n-gram analysis.

The analysis conducted in this article resulted in a predictor of amyloidogenicity AmyloGram, available as a *R* package (https://cran.r-project.org/package=AmyloGram) and a web server (www.smorfland.uni.wroc.pl/shiny/AmyloGram/).

## How to reproduce the main part of the analysis?

Source **analysis.R**. Be warned that computations are time consuming.

## How to generate results and plots for publication?

Run **analysis_publication.R**. This step can be executed even without performing the full analysis, because the partial results are avaible in the repository.

## Repository structure

### benchmark

Raw data and results of the benchmark.

### data

Processed data used in the study: 
* amyloid and non-amyloid sequences derived from [AmyLoad](comprec-lin.iiar.pwr.edu.pl/amyload) database on 09.01.2016.
* amyloid and non-amyloid benchmark data from *pep424* dataset.
* *tableA.csv* normalized table A from supplemental materials of [Wozniak2014](https://dx.doi.org/10.1007/s00894-014-2497-9).
* *AA_index_mk2.csv* amino acid properties chosen by MK.

### data exploration

Scripts used in the exploratory phase of data analysis. Usually undocumented. All exploratory scripts that were found useful were later incorporated into the analysis and belog to the **functions** directory. 

### functions

All functions necessary to repeat the analysis. Functions are documented using roxygen standard. 

### varia

Additional source files as the logo of the AmyloGram package.

## R session information

All scripts used in this study are compatible with following versions of R and necessary packages:

**R version 3.4.1 (2017-06-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=pl_PL.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=pl_PL.UTF-8_, _LC_COLLATE=pl_PL.UTF-8_, _LC_MONETARY=pl_PL.UTF-8_, _LC_MESSAGES=pl_PL.UTF-8_, _LC_PAPER=pl_PL.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=pl_PL.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_grid_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_bindrcpp(v.0.2)_, _dplyr(v.0.7.1)_, _reshape2(v.1.4.2)_, _xtable(v.1.8-2)_, _seqinr(v.3.3-6)_, _biogram(v.1.5)_, _slam(v.0.1-40)_, _ggrepel(v.0.6.5)_, _gtable(v.0.2.0)_, _gridExtra(v.2.2.1)_, _ggplot2(v.2.2.1)_ and _pander(v.0.6.0)_

**loaded via a namespace (and not attached):** 
_Rcpp(v.0.12.11)_, _bindr(v.0.1)_, _magrittr(v.1.5)_, _partitions(v.1.9-18)_, _munsell(v.0.4.3)_, _bit(v.1.1-12)_, _colorspace(v.1.3-2)_, _R6(v.2.2.2)_, _rlang(v.0.1.1)_, _stringr(v.1.2.0)_, _plyr(v.1.8.4)_, _tools(v.3.4.1)_, _assertthat(v.0.2.0)_, _ade4(v.1.7-6)_, _lazyeval(v.0.2.0)_, _digest(v.0.6.12)_, _tibble(v.1.3.3)_, _gmp(v.0.5-13.1)_, _glue(v.1.1.1)_, _stringi(v.1.1.5)_, _entropy(v.1.2.1)_, _compiler(v.3.4.1)_, _scales(v.0.4.1)_, _combinat(v.0.0-8)_, _polynom(v.1.3-9)_ and _pkgconfig(v.2.0.1)_