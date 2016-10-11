# How to reproduce the main part of the analysis?

Source **analysis.R**. Be warned that computations are time consuming.

# How to generate results and plots for publication?

Run **analysis_publication.R**. This step can be executed even without performing the full analysis, because the partial results are avaible in the repository.

# Repository structure

## benchmark

Raw data and results of the benchmark.

## data

Processed data used in the study: 
* amyloid and non-amyloid sequences derived from [AmyLoad](comprec-lin.iiar.pwr.edu.pl/amyload) database on 09.01.2016.
* amyloid and non-amyloid benchmark data from *pep424* dataset.
* *tableA.csv* normalized table A from supplemental materials of [Wozniak2014](https://dx.doi.org/10.1007/s00894-014-2497-9).
* *AA_index_mk2.csv* amino acid properties chosen by MK.

## data exploration

Scripts used in the exploratory phase of data analysis. Usually undocumented. All exploratory scripts that were found useful were later incorporated into the analysis and belog to the **functions** directory. 

## functions

All functions necessary to repeat the analysis. Functions are documented using roxygen standard. 

## varia

Additional source files as the logo of the AmyloGram package.

# R session information

All scripts used in this study are compatible with following versions of R and necessary packages:

**R version 3.3.1 (2016-06-21)**

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

* ranger(v.0.5.0) 
* dplyr(v.0.5.0) 
* reshape2(v.1.4.1) 
* xtable(v.1.8-2) 
* seqinr(v.3.3-1) 
* biogram(v.1.3) 
* slam(v.0.1-38) 
* gridExtra(v.2.2.1) 
* ggplot2(v.2.1.0) 


**loaded via a namespace (and not attached):** 

* Rcpp(v.0.12.7) 
* magrittr(v.1.5) 
* partitions(v.1.9-18) 
* bit(v.1.1-12) 
* munsell(v.0.4.3) 
* colorspace(v.1.2-6) 
* R6(v.2.2.0) 
* stringr(v.1.1.0) 
* plyr(v.1.8.4) 
* tools(v.3.3.1) 
* gtable(v.0.2.0) 
* DBI(v.0.5-1) 
* htmltools(v.0.3.5) 
* digest(v.0.6.10) 
* ade4(v.1.7-4) 
* assertthat(v.0.1) 
* tibble(v.1.2) 
* gmp(v.0.5-12) 
* evaluate(v.0.9) 
* rmarkdown(v.1.0) 
* entropy(v.1.2.1) 
* stringi(v.1.1.2) 
* scales(v.0.4.0) 
* combinat(v.0.0-8) 
* polynom(v.1.3-8) 