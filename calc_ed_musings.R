a <- structure(list(`1` = "p", `2` = c("f", "i", "w", "y"), 
                    `3` = c("a", "c", "d", "e", "g", "h", "k", "l", "m", "n", "q", "r", "s", "t", "v")), 
               .Names = c("1", "2", "3"))

b <- structure(list(`1` = c("f", "r", "w", "y"), `2` = c("c", "i", "l", "t", "v"), 
                    `3` = c("a", "d", "e", "g", "h", "k", "m", "n", "p", "q", "s")), 
               .Names = c("1", "2", "3"))

calc_ed(a, b)

# preserve for later - code for recovering 
# traits <- create_traits_combination(ftraits)
# 
# traits_table <- do.call(rbind, lapply(traits, function(i) 
#   cbind(i, matrix(NA, ncol = 17 - ncol(i), nrow = nrow(i)))))
# 
# id_length_table <- names(aa_groups) %>% 
#   sapply(function(i) substr(i, 3, nchar(i))) %>%
#   strsplit("K") %>%
#   do.call(rbind, .) %>%
#   apply(2, as.numeric) 
# 
# 
# aa_nprop <- normalize_properties()
# 
# a_prop <- aa_nprop[na.omit(traits_table[50, ]), , drop = FALSE]

a_prop <- structure(c(0.112469437652812, 0.405737704918033, 0.711491442542787, 
                      0.332991803278688, 0.327628361858191, 0.156762295081967, 0.256723716381418, 
                      0.00614754098360652, 0.312958435207824, 0.905737704918033, 0.440097799511002, 
                      0.145491803278689, 0.36919315403423, 0.055327868852459, 0, 0.262295081967213, 
                      0.562347188264059, 0.559426229508197, 0.454767726161369, 1, 0.454767726161369, 
                      0.941598360655738, 0.535452322738386, 0, 0.540342298288509, 0.787909836065574, 
                      0.709046454767726, 0.968237704918033, 0.320293398533007, 0.117827868852459, 
                      0.15158924205379, 0.137295081967213, 0.264058679706601, 0.305327868852459, 
                      1, 0.961065573770492, 0.728606356968215, 0.648565573770492, 0.342298288508558, 
                      0.88422131147541), .Dim = c(2L, 20L), 
                    .Dimnames = list(c("CHAM820101", "NISK860101"), 
                                     c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")))

colnames(a_prop) <- tolower(colnames(a_prop))

#b_prop <- aa_nprop[na.omit(traits_table[ao, ]), , drop = FALSE]

# must have unified lists of features

coords_a <- lapply(a, function(single_subgroup) rowMeans(a_prop[, single_subgroup, drop = FALSE]))
coords_b <- lapply(b, function(single_subgroup) rowMeans(a_prop[, single_subgroup, drop = FALSE]))

# finished - distance as the minimum difference in properties between groups
# how to deal with encodings with the unequal number of groups
sapply(coords_a, function(single_coords_a) {
  distances <- sapply(coords_b, function(single_coords_b) 
    #vector of distances between groups
    sqrt(sum((single_coords_a - single_coords_b)^2))
  )
  #c(dist = min(distances), id = unname(which.min(distances)))
  min(distances)
})

dat_a <- data.frame(enc = "a", do.call(rbind, coords_a))
dat_b <- data.frame(enc = "b", do.call(rbind, coords_b))

dat <- data.frame(do.call(rbind, lapply(1L:nrow(dat_a), function(id) 
  data.frame(id = id, rbind(do.call(rbind, lapply(1L:3, function(dummy) 
    dat_a[id, , drop = FALSE])),
    dat_b)))), pair = c(paste0("d", 1L:3), paste0("d", 1L:3)))

colnames(dat) <- c("id", "enc", "f1", "f2", "pair")
dat[["id"]] <- paste0("Encoding A, subgroup ", dat[["id"]])


library(ggplot2)
ggplot(dat, aes(x = f1, y = f2, colour = pair)) +
  geom_line() +
  geom_point(aes(x = f1, y = f2, colour = enc), size = 10) + 
  facet_wrap(~ id) + 
  scale_color_brewer(palette="Dark2", guide = "none")



my_group <- list(`1` = c("c", "f", "h", "i", "l", "v", "w"),
                 `2` = c("a", "d", "g", "k", "n", "p", "q", "r"),
                 `3` = c("e", "m", "s", "t", "y"))

all_dists <- pbsapply(aa_groups, function(i) calc_ed(i, my_group))

tmp <- aa_groups[all_dists == 1][[1]]

calc_ed(tmp, my_group)
