b <- aa_groups[[5]]
a <- aa_groups[[5000]]

library(biogram)
calc_ed(a, b)

traits <- create_traits_combination(vtraits)

traits_table <- do.call(rbind, lapply(traits, function(i) 
  cbind(i, matrix(NA, ncol = 17 - ncol(i), nrow = nrow(i)))))

id_length_table <- names(aa_groups) %>% 
  sapply(function(i) substr(i, 3, nchar(i))) %>%
  strsplit("K") %>%
  do.call(rbind, .) %>%
  apply(2, as.numeric) 


traits_table[5000, ]

aa_nprop <- normalize_properties()

a_prop <- aa_nprop[na.omit(traits_table[5000, ]), , drop = FALSE]
b_prop <- aa_nprop[na.omit(traits_table[5, ]), , drop = FALSE]
