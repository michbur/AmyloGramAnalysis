source("choose_properties.R")

library(dplyr)
library(reshape2)

plot_values <- prop_MK %>% 
  select(X) %>% 
  unlist %>% 
  c(., 545L:550) %>% 
  slice(data.frame(aa_nprop), .) %>%
  t 

mplot_values <- melt(plot_values)

colnames(mplot_values) <- c("aa", "id", "value")
mplot_values[["id"]] <- as.factor(mplot_values[["id"]])


plot_names <- prop_MK %>% select(name) %>% unlist %>% as.character %>% c(., add_names)
plot_id <- 1L:length(plot_names)
names(plot_id) <- plot_names


save(plot_values, mplot_values, plot_id, file = "./property_app/properties.RData")

#input$checkProp
# tmp <- c(1, 2, 3)
# 
# library(ggdendro)
# cl <- hclust(dist(plot_values[, as.numeric(tmp)]))
# gr <- cutree(cl, k = 3)
# 
# ddata <- dendro_data(dhc, type = "rectangle")
# p <- ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   coord_flip() + 
#   scale_y_reverse(expand = c(0.2, 0))

#PCA
tmp <- c(17, 18, 19, 20)
plot_values[, as.numeric(tmp)]
pca_res <- princomp(plot_values[, as.numeric(tmp)])
pca_res[["loadings"]]

ggplot(data.frame(pca_res[["scores"]]), aes(x = Comp.1, y = Comp.2, colour = Comp.3, label = rownames(pca_res[["scores"]]))) +
  geom_point(size = 5) +
  scale_colour_gradient2(low = "red", high = "yellow", mid = "orange") +
  geom_text(size = 5)

#correlation plot


tmp <- 1L:26
nms <- names(plot_id)[as.numeric(tmp)]
br_nms <- sapply(nms, function(single_nm) {
  len <- nchar(single_nm)
  space_pos <-  which(strsplit(single_nm, "")[[1]] == " ")
  break_pos <- space_pos[which.min(abs(len/2 - space_pos))]
  paste0(substr(single_nm, 1, break_pos - 1), "\n", substr(single_nm, break_pos + 1, len))
})

corm <- melt(cor(plot_values[, as.numeric(tmp)]))
corm[["Var1"]] <- factor(corm[["Var1"]], labels = br_nms)
corm[["Var2"]] <- factor(corm[["Var2"]], labels = br_nms)

corm[["value_fac"]] <- round(corm[["value"]], 4)

ggplot(data = corm, aes(x=Var1, y=Var2, fill=value, label = value_fac)) + 
  geom_tile() +
  geom_text()

corm <- melt(cor(plot_values[, as.numeric(tmp)]))
nms <- names(plot_id)[as.numeric(tmp)]
corm[["Var1"]] <- factor(corm[["Var1"]], labels = nms)
corm[["Var2"]] <- factor(corm[["Var2"]], labels = nms)

lapply(levels(corm[["Var1"]]), function(i)
  corm[corm[["Var1"]] == i & corm[["Var2"]] != i & abs(corm[["value"]]) > 0.8, ])
