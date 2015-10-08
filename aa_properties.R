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

#properties below are hydrophobicity scale that should be reversed
aa_nprop[c(252, 519, 543, 544), ] <- 1 - aa_nprop[c(252, 519, 543, 544), ]

#key for selecting properties - new (younger than XX years properties)

traits <- list(size = 515,
               hydroph = c(494, 529, 528),
               solvent = c(319, 211, 512),
               polarity = c(22, 321),
               interactivity = 545:550)