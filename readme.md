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
