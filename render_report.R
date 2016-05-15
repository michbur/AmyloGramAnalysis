library(rmarkdown)
library(knitr)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dokumenty/gits/amyloid_ngram2/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "D:/GitHub/amyloid_ngram2"

opts_knit$set(root.dir = pathway)

render("./reports/methods.Rmd")

render("./reports/reg33_analysis/reg33.Rmd")
