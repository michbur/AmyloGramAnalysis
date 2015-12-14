library(rmarkdown)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dokumenty/gits/amyloid_ngram2/"

opts_knit$set(root.dir = pathway)

render("./reports/report1.Rmd")

