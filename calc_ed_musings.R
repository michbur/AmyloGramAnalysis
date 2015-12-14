bi <- 50
ai <- 5000
b <- aa_groups[[bi]]
a <- aa_groups[[ai]]

library(biogram)
calc_ed(a, b)

traits <- create_traits_combination(ftraits)

traits_table <- do.call(rbind, lapply(traits, function(i) 
  cbind(i, matrix(NA, ncol = 17 - ncol(i), nrow = nrow(i)))))

id_length_table <- names(aa_groups) %>% 
  sapply(function(i) substr(i, 3, nchar(i))) %>%
  strsplit("K") %>%
  do.call(rbind, .) %>%
  apply(2, as.numeric) 


traits_table[5000, ]

aa_nprop <- normalize_properties()

a_prop <- aa_nprop[na.omit(traits_table[ai, ]), , drop = FALSE]
colnames(a_prop) <- tolower(colnames(a_prop))

#b_prop <- aa_nprop[na.omit(traits_table[ao, ]), , drop = FALSE]

# must have unified lists of features

coords_a <- lapply(a, function(single_subgroup) rowMeans(a_prop[, single_subgroup, drop = FALSE]))
coords_b <- lapply(b, function(single_subgroup) rowMeans(a_prop[, single_subgroup, drop = FALSE]))

# finished - distance as the minimum difference in properties between groups
# how to deal with encodings with the unequal number of groups
sum(sapply(coords_a, function(single_coords_a) {
  distances <- sapply(coords_b, function(single_coords_b) 
    #vector of distances between groups
    sqrt(sum((single_coords_a - single_coords_b)^2))
  )
  #c(dist = min(distances), id = unname(which.min(distances)))
  min(distances)
}))
