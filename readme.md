# Read me

This repository contains data and code necessary to reproduce analysis from the article: [Burdukiewicz M, Sobczyk P, Rödiger S, Duda-Madej A, Mackiewicz P, Kotulska M. (2017) Amyloidogenic motifs revealed by n-gram analysis (Scientific Reports 7 (1), 12961)](https://www.nature.com/articles/s41598-017-13210-9).

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

**attached base packages:** 

* grid 
* stats 
* graphics 
* grDevices 
* utils 
* datasets 
* methods 
* base 


**other attached packages:** 

* bindrcpp(v.0.2) 
* dplyr(v.0.7.1) 
* reshape2(v.1.4.2) 
* xtable(v.1.8-2) 
* seqinr(v.3.3-6) 
* biogram(v.1.5) 
* slam(v.0.1-40) 
* ggrepel(v.0.6.5) 
* gtable(v.0.2.0) 
* gridExtra(v.2.2.1) 
* ggplot2(v.2.2.1) 


**loaded via a namespace (and not attached):** 

* Rcpp(v.0.12.11) 
* bindr(v.0.1) 
* magrittr(v.1.5) 
* partitions(v.1.9-18) 
* bit(v.1.1-12) 
* munsell(v.0.4.3) 
* colorspace(v.1.3-2) 
* R6(v.2.2.2) 
* rlang(v.0.1.1) 
* stringr(v.1.2.0) 
* plyr(v.1.8.4) 
* tools(v.3.4.1) 
* digest(v.0.6.12) 
* assertthat(v.0.2.0) 
* ade4(v.1.7-6) 
* lazyeval(v.0.2.0) 
* tibble(v.1.3.3) 
* gmp(v.0.5-13.1) 
* glue(v.1.1.1) 
* entropy(v.1.2.1) 
* stringi(v.1.1.5) 
* pander(v.0.6.0) 
* compiler(v.3.4.1) 
* scales(v.0.4.1) 
* combinat(v.0.0-8) 
* polynom(v.1.3-9) 
* pkgconfig(v.2.0.1) 