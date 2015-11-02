library(seqinr)
library(dplyr)
library(pbapply)

paste_enc <- function(x)
  paste0(sapply(x, paste0, collapse = ""), collapse = "_")

source("choose_properties.R")

grouping_properties <- aa_nprop[unlist(ftraits), ]

#vector of traits
vtraits <- unlist(ftraits)

#all combinations of traits
all_traits_combn_list <- pblapply(1L:length(vtraits), function(i)
  t(combn(vtraits, i)))
# number of combinations: sum(sapply(all_traits_combn_list, nrow))

#create encodings
all_aa_groups <- pblapply(3L:6, function(single_k) {
  res <- unlist(lapply(all_traits_combn_list, function(all_traits_combn)
    vapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
      cl <- t(aa_nprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE]) %>%
        dist %>%
        hclust(method = "ward.D2")
      #cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE])))
      gr <- cutree(cl, k = single_k)
      names(gr) <- tolower(names(gr))
      agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
      #inside encodings, amino acids are ordered alphabetically
      agg_gr <- lapply(agg_gr, sort)
      #groups are sorted by their length
      paste_enc(agg_gr[order(lengths(agg_gr))])
    }, "a")))
  names(res) <- paste0("ID", 1L:length(res), "K", single_k)
  res
})

#get indices of unique encodings
aa_id <- lapply(all_aa_groups, function(i) !duplicated(i))

#remove from aa_groups redundant encodings
aa_groups <- unlist(lapply(1L:length(aa_id), function(i) {
  all_aa_groups[[i]][aa_id[[i]]]
}), recursive = FALSE)

#add as a benchmark two encodings from the literature
aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))

aa_groups <- c(aa1 = paste_enc(aa1), aa2 = paste_enc(aa2), aa_groups)

save(aa_groups, file = "aa_groups.RData")
