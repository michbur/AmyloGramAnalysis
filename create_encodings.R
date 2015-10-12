library(seqinr)
library(dplyr)

source("choose_properties.R")

grouping_properties <- aa_nprop[unlist(traits), ]

#vector of traits
vtraits <- unlist(traits)

#all combinations of traits
all_traits_combn_list <- lapply(1L:15, function(i)
  t(combn(vtraits, i)))

#table of features-encodings (features on which encodings are based)
fenc_table <- do.call(rbind, lapply(all_traits_combn_list, function(i)
  cbind(i, matrix(NA, nrow = nrow(i), ncol = length(vtraits) - ncol(i)))))
rownames(fenc_table) <- paste0("ID", 1L:nrow(enc_table))

#create encodings
aa_groups <- lapply(3L:6, function(single_k) {
  res <- lapply(all_traits_combn_list, function(all_traits_combn)
    lapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
      cl <- t(aa_nprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE]) %>%
        dist %>%
        hclust(method = "ward.D2")
      #cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE])))
      gr <- cutree(cl, k = single_k)
      names(gr) <- tolower(names(gr))
      agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
      names(agg_gr) <- 1L:length(agg_gr)
      agg_gr
    })) %>% unlist(recursive = FALSE) 
  names(res) <- paste0("ID", 1L:length(res), "K", single_k)
  res
})

#get indices of unique encodings
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

aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))

aa_groups <- c(list(aa2), list(aa1), aa_groups)