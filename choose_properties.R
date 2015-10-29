#' @description Normalize amino acid properties, add properties from Wozniak 2014, 
#' create a list of traits for encodings

library(seqinr)
library(dplyr)

data("aaindex")

aa_props <- sapply(aaindex, function(i) i[["I"]])
tableA <- read.table("tableA.csv", sep = ";", dec = ".", head = TRUE)

#names for additional properties
add_names <- c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")


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

prop_MK <- cbind(prop_MK, years = years) %>% filter(years >= 1980)

traits <- c(prop_MK[["X"]], 545:550)

traits_names <- c(sapply(prop_MK[["X"]], function(i) aaindex[[i]][["D"]]), 
                  add_names)
#final traits
ftraits <- traits[-c(6, 13, 18, 19, 20, 21, 22, 24, 25)]
