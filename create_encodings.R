library(seqinr)
library(dplyr)

source("choose_properties.R")

grouping_properties <- aa_nprop[unlist(traits), ]

all_traits_combn_list <- list(expand.grid(traits))

aa_groups <- unlist(lapply(all_traits_combn_list, function(all_traits_combn)
  lapply(3L:6, function(single_k)
    lapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
      cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, ]), ])))
      gr <- cutree(cl, k = single_k)
      names(gr) <- tolower(names(gr))
      agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
      names(agg_gr) <- 1L:length(agg_gr)
      agg_gr
    }))), recursive = FALSE)

aa_id <- lapply(aa_groups, function(j) {
  sort_gr <- lapply(j, function(i) {
    res <- sapply(i, sort)
    res[order(lengths(res))]
  })
  !duplicated(t(sapply(sort_gr, unlist)))
})

aa_groups <- unlist(lapply(1L:length(aa_id), function(i) {
  aa_groups[[i]][aa_id[[i]]]
}), recursive = FALSE)

aa_groups <- sapply(aa_groups, function(i) {
  res <- sapply(i, sort)
  res[order(lengths(res))]
})

groups_summary <- cbind(chosen = !duplicated(t(sapply(aa_groups, unlist))), do.call(rbind, lapply(3L:6, function(i) 
  cbind(n = rep(i, nrow(all_traits_combn_list[[1]])), all_traits_combn_list[[1]]))))

aa_groups <- aa_groups[!duplicated(t(sapply(aa_groups, unlist)))]

aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))


aa_groups <- c(list(aa2), list(aa1), aa_groups)