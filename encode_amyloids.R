load("aa_groups.RData")

#recover encoding from the string
recover_enc <- function(x)
  lapply(strsplit(x,  split = "_"), strsplit, split = "")
