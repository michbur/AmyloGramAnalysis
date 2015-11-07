load("aa_groups.RData")

#recover encoding from the string
vector2list <- function(x){
  require(dplyr)
  lapply(x, function(i) {
    pasted_group <- strsplit(i, "_", fixed = TRUE)
    lapply(pasted_group, function(j) strsplit(j, ""))
  }) %>% lapply(first)
}
