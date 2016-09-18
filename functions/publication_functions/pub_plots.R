# plot data ----------------------------------------------

amyloids_plot <- select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
  mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
  select(AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, et, enc_adj) %>%
  rbind(select(full_alphabet, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>% 
          mutate(et = "full alphabet", enc_adj = 0)) %>%
  mutate(len_range = factor(len_range, levels = c("[5,6]", "(6,10]", "(10,15]", "(15,25]")),
         pos = factor(pos, labels = paste0("Training peptide length: ", 
                                           c("6", "6-10", "6-15"))),
         et = factor(et, labels = c("Encoding", "Best-performing encoding",
                                    "Standard encoding", "Full alphabet"))) %>%
  mutate(len_range = factor(len_range, 
                            labels = paste0("Test peptide length: ", c("6 ", "7-10", "11-15", "16-25"))),
         et2 = ifelse(enc_adj == 1L, "Standard encoding (Kosiol et al., 2004)", as.character(et)),
         et2 = ifelse(enc_adj == 2L, "Standard encoding (Melo and Marti-Renom, 2006)", as.character(et2)),
         et2 = factor(et2, levels = c("Encoding", 
                                      "Best-performing encoding", 
                                      "Full alphabet", 
                                      "Standard encoding (Kosiol et al., 2004)", 
                                      "Standard encoding (Melo and Marti-Renom, 2006)")),
         et = et2)

# write.csv(amyloids_plot, file = "./results/amyloid_plot_data.csv")

# Fig 1 all encodings sens/spec  ----------------------------------------

sesp_dat <- amyloids_plot
levels(sesp_dat[["pos"]]) <- c("Training peptide\nlength: 6", "Training peptide\nlength: 6-10", 
                               "Training peptide\nlength: 6-15")

# sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et)) +
#   #geom_density_2d(color = "grey", contour = TRUE) +
#   #stat_density2d(aes(fill=..level..), color = "grey", contour = TRUE, geom="polygon") +
#   scale_alpha_continuous(range = c(0.35, 1)) +
#   scale_y_continuous("Mean sensitivity") +
#   scale_x_continuous("Mean specificity") +
#   geom_point(data = sesp_dat, aes(shape = et)) +
#   geom_point(data = filter(sesp_dat, et != "Encoding"), aes(shape = et)) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
#   scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
#   scale_color_manual("", values = c("grey", "red", "green", "blue", "blue"), drop = FALSE) +
#   scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) +
#   facet_grid(pos ~ len_range) +
#   my_theme 


sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et)) +
  geom_bin2d(bins = 30, color = "black") + 
  scale_fill_continuous("Count", low = "gray83", high = "gray44") +
  scale_y_continuous("Mean sensitivity") +
  scale_x_continuous("Mean specificity") +
  geom_point(data = droplevels(filter(sesp_dat, et != "Encoding")),
             aes(x = Spec_mean, y = Sens_mean, color = et2, shape = et2), size = 4) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), 
         fill = guide_colorbar(barwidth = unit(50, "line"))) +
  scale_shape_manual("", values = c(16, 15, 17, 17), drop = FALSE) +
  scale_color_manual("", values = c("red", "green", "blue", "blue"), drop = FALSE) +
  scale_size_manual("", values = c(1, 1, 1.5, 1.5), drop = FALSE) +
  facet_grid(pos ~ len_range) +
  my_theme


png("./publication/figures/sesp_plot.png", height = 4, width = 6.5, unit = "in", res = 200)
#cairo_ps("./pub_figures/sesp_plot.eps", height = 4, width = 8)
# should be eps, but it's too big for overleaf
print(sesp_plot)
dev.off()

# Fig 2 AUC boxplot  ----------------------------------------

AUC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = AUC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1, outlier.size = 1) +
  geom_point(data = filter(amyloids_plot, et2 != "Encoding"), 
             aes(x = len_range, y = AUC_mean, color = et2, shape = et2, size = et2)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean AUC") +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "red", "green", "blue", "blue"), drop = FALSE) +
  scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip() 

cairo_ps("./publication/figures/AUC_boxplot.eps", height = 3.5, width = 6.5)
#png("./pub_figures/AUC_boxplot.png", height = 648, width = 648)
print(AUC_boxplot)
dev.off()

# Fig 3 MCC boxplot  ----------------------------------------

MCC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = MCC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1) +
  geom_point(data = filter(amyloids_plot, et != "Encoding"), 
             aes(x = len_range, y = MCC_mean, color = et, shape = et)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean MCC") +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "red", "blue", "green"), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip()

# cairo_ps("./publication/figures/MCC_boxplot.eps", height = 9, width = 8)
# #png("./pub_figures/MCC_boxplot.png", height = 648, width = 648)
# MCC_boxplot
# dev.off()

# Fig 4 properties  ----------------------------------------

ggplot(best_enc_props, aes(x = as.factor(id), y = value, label = aa)) +
  geom_text(position = "dodge") +
  facet_wrap(~ gr, ncol = 2)


# Fig 5 n-grams  ----------------------------------------

gr_aa <- group_by(best_enc_aa, id) %>% 
  summarise(gr = paste0("(", paste0(aa, collapse = ", "), ")"))

ngram_freq_plot <- mutate(ngram_freq, decoded_name = gsub("_", " - ", decoded_name)) %>%
  mutate(decoded_name = factor(decoded_name, levels = as.character(decoded_name)),
         amyloid = diff_freq > 0) %>%
  melt() %>%
  filter(variable %in% c("pos", "neg")) %>%
  droplevels %>%
  mutate(variable = factor(variable, labels = c("Amyloid", "Non-amyloid")))

for(i in 1L:6)
  levels(ngram_freq_plot[["decoded_name"]]) <- gsub(as.character(i), gr_aa[i, "gr"],
    levels(ngram_freq_plot[["decoded_name"]]))
                                                    

ngram_plot <- ggplot(ngram_freq_plot, aes(x = decoded_name, y = value)) +
  geom_bar(aes(fill = variable), position = "dodge", stat = "identity") +
  geom_point(data = group_by(ngram_freq_plot, decoded_name)  %>% filter(value == max(value)),
             aes(y = value + 0.002, shape = association)) +
  scale_fill_manual("", values = c("gold", "darkmagenta")) +
  scale_shape_manual("Motif:", breaks = c("Amyloidogenic", "Non-amyloidogenic"), values = c(16, 17, NA)) +
  scale_y_continuous("Frequency") +
  scale_x_discrete("") +
  coord_flip() +
  my_theme 

# in case we need to get n-grams in a tabular format
#writeLines(as.character(ngram_freq_plot[["decoded_name"]]), "n_gramy_Ania.txt")

cairo_ps("./publication/figures/ngrams.eps", height = 8, width = 2.5)
print(ngram_plot)
dev.off()

# Fig 6 encoding distance  ----------------------------------------

# careful - check if similarity index is used instead of the encoding distance

si_dat <- si_dat %>% 
  mutate(et2 = ifelse(enc_adj == 1L, "Standard encoding (Kosiol et al., 2004)", as.character(et)),
         et2 = ifelse(enc_adj == 2L, "Standard encoding (Melo and Marti-Renom, 2006)", as.character(et2)),
         et2 = factor(et2, levels = c("Encoding", 
                                      "Best-performing encoding", 
                                      "Full alphabet", 
                                      "Standard encoding (Kosiol et al., 2004)", 
                                      "Standard encoding (Melo and Marti-Renom, 2006)")),
         et = et2)

write.csv2(si_dat, row.names = FALSE, file = "./results/si_dat.csv")

# ed_AUC_plot <- ggplot(si_dat, aes(x=si, y=AUC_mean, color=et, shape = et)) + 
#   geom_point() +
#   scale_color_manual("", values = c("grey", "red", "blue", "green")) +
#   scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
#   xlab("Similarity") +
#   ylab("AUC") +
#   my_theme +
#   geom_point(data = filter(si_dat, et != "Encoding"), 
#              aes(x = si, y = AUC_mean, color = et)) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 5)) +
#   scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
#   scale_color_manual("", values = c("grey", "red", "green", "blue", "blue"), drop = FALSE) +
#   scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) 

# Fig 6 alternative (similarity index)  ----------------------------------------

si_dat <- si_dat %>% 
  mutate(et2 = ifelse(enc_adj == 1L, "Standard encoding (Kosiol et al., 2004)", as.character(et)),
         et2 = ifelse(enc_adj == 2L, "Standard encoding (Melo and Marti-Renom, 2006)", as.character(et2)),
         et2 = factor(et2, levels = c("Encoding", 
                                      "Best-performing encoding", 
                                      "Full alphabet", 
                                      "Standard encoding (Kosiol et al., 2004)", 
                                      "Standard encoding (Melo and Marti-Renom, 2006)")),
         et = et2)

# si_AUC_plot <- ggplot(si_dat, aes(x=si, y=AUC_mean, color=et, shape = et)) + 
#   geom_point() +
#   scale_color_manual("", values = c("grey", "red", "blue", "green")) +
#   scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
#   xlab("Similarity index") +
#   ylab("AUC") +
#   my_theme +
#   geom_point(data = filter(si_dat, et != "Encoding"), 
#              aes(x = si, y = AUC_mean, color = et)) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 5)) +
#   scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
#   scale_color_manual("", values = c("grey", "red", "green", "blue", "blue"), drop = FALSE) +
#   scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) 

si_AUC_plot <- ggplot(si_dat, aes(x=si, y=AUC_mean)) + 
  geom_bin2d(bins = 30, color = "black") + 
  scale_fill_continuous("Count", low = "gray83", high = "gray44") +
  xlab("Similarity") +
  ylab("AUC") +
  my_theme +
  geom_point(data = droplevels(filter(si_dat, et != "Encoding")),
             aes(x = si, y = AUC_mean, color = et2, shape = et2), size = 4) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), 
         fill = guide_colorbar(barwidth = unit(50, "line"))) +
  scale_shape_manual("", values = c(16, 15, 17, 17), drop = FALSE) +
  scale_color_manual("", values = c("red", "green", "blue", "blue"), drop = FALSE) +
  scale_size_manual("", values = c(1, 1, 1.5, 1.5), drop = FALSE)

cairo_ps("./publication/figures/ed_AUC.eps", height = 4, width = 3)
print(si_AUC_plot)
dev.off()


# save(amyloids_plot, best_enc_props, ngram_freq_plot, ed_dat,
#      file = "./presentation/presentation.RData")
