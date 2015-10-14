#' @description Normalize amino acid properties, add properties from Wozniak 2014, 
#' create a list of traits for encodings

library(seqinr)
library(dplyr)

data("aaindex")

aa_props <- sapply(aaindex, function(i) i[["I"]])
tableA <- read.table("tableA.csv", sep = ";", dec = ".", head = TRUE)

aa_nprop <- t(apply(cbind(aa_props, tableA[tableA[["X"]] %>% as.character %>% aaa %>% order, 2L:7]), 2, function(i) {
  res <- i - min(i, na.rm = TRUE)
  res/max(res, na.rm = TRUE)
}))

colnames(aa_nprop) <- a(colnames(aa_nprop))

#properties below are hydrophobicity scales that should be reversed
aa_nprop[c(252, 519, 543, 544), ] <- 1 - aa_nprop[c(252, 519, 543, 544), ]

#key for selecting properties - new (younger than 1980 years properties)
prop_MK <- read.csv2("AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1]) %>%
  as.numeric

prop_MK <- cbind(prop_MK, years = years) %>% filter(years > 1980)

traits <- c(prop_MK[["X"]], 545:550)
