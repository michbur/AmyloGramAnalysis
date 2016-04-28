# plot data ----------------------------------------------

amyloids_plot <- select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
  mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
  select(AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, et, enc_adj) %>%
  rbind(select(full_alphabet, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>% 
          mutate(et = "full alphabet", enc_adj = 0)) %>%
  mutate(len_range = factor(len_range, levels = c("[5,6]", "(6,10]", "(10,15]", "(15,25]")),
         pos = factor(pos, labels = paste0("Training peptide length: ", 
                                           c("6", "6-10", "6-15"))),
         et = factor(et, labels = c("Reduced alphabet", "Best performing reduced alphabet",
                                    "Reduced alphabet from literature", "Full alphabet"))) %>%
  mutate(len_range = factor(len_range, 
                            labels = paste0("Testing peptide length: ", c("6 ", "7-10", "11-15", "16-25"))))

# Fig 1 all encodings sens/spec  ----------------------------------------

sesp_dat <- amyloids_plot
levels(sesp_dat[["pos"]]) <- c("Training peptide\nlength: 6", "Training peptide\nlength: 6-10", 
                               "Training peptide\nlength: 6-15")

sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et, shape = et)) +
  geom_point() +
  scale_color_manual("", values = c("grey", "red", "blue", "green")) +
  scale_shape_manual("", values = c(1, 16, 15, 15)) +
  scale_y_continuous("Mean sensitivity") +
  scale_x_continuous("Mean specificity") +
  facet_grid(pos ~ len_range) +
  my_theme +
  geom_point(data = filter(sesp_dat, et != "Reduced alphabet"), 
             aes(x = Spec_mean, y = Sens_mean, color = et))

png("./publication/figures/sesp_plot.png", height = 4, width = 6.5, unit = "in", res = 200)
#cairo_ps("./pub_figures/sesp_plot.eps", height = 4, width = 8)
# should be eps, but it's too big for overleaf
print(sesp_plot)
dev.off()

# Fig 2 AUC boxplot  ----------------------------------------

AUC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = AUC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1) +
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
             aes(x = len_range, y = AUC_mean, color = et, shape = et)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean AUC") +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "red", "blue", "green"), drop = FALSE) +
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
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
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

ngram_freq_plot <- mutate(ngram_freq, decoded_name = gsub("_", "|", decoded_name)) %>%
  mutate(decoded_name = factor(decoded_name, levels = as.character(decoded_name)),
         amyloid = diff_freq > 0) %>%
  melt() %>%
  filter(variable %in% c("pos", "neg")) %>%
  droplevels %>%
  mutate(variable = factor(variable, labels = c("Amyloid", "Non-amyloid")))

ngram_plot <- ggplot(ngram_freq_plot, aes(x = decoded_name, y = value, fill = variable)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual("", values = c("gold", "darkmagenta")) +
  scale_y_continuous("Frequency", limits = c(0, 0.95)) +
  scale_x_discrete("") +
  coord_flip() +
  my_theme + 
  theme(panel.grid.major.y = element_line(color = "NA")) 

cairo_ps("./publication/figures/ngrams.eps", height = 8.5, width = 2.5)
print(ngram_plot)
dev.off()

# Fig 6 encoding distance  ----------------------------------------

ed_AUC_plot <- ggplot(ed_dat, aes(x=ed, y=AUC_mean, color=et, shape = et)) + 
  geom_point() +
  scale_color_manual("", values = c("grey", "red", "blue", "green")) +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  xlab("Normalized encoding distance") +
  ylab("AUC") +
  my_theme +
  geom_point(data = filter(ed_dat, et != "Reduced alphabet"), 
             aes(x = ed, y = AUC_mean, color = et)) +
  guides(color=guide_legend(ncol=2))

cairo_ps("./publication/figures/ed_AUC.eps", height = 4, width = 3)
print(ed_AUC_plot)
dev.off()

# save(amyloids_plot, best_enc_props, ngram_freq_plot, ed_dat,
#      file = "./presentation/presentation.RData")
