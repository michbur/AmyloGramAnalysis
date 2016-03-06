# plot data ----------------------------------------------

amyloids_plot <- select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
  mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
  select(AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, et, enc_adj) %>%
  rbind(select(full_alphabet, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>% 
          mutate(et = "full alphabet", enc_adj = 0)) %>%
  mutate(len_range = factor(len_range, levels = c("[5,6]", "(6,10]", "(10,15]", "(15,25]")),
         pos = factor(pos, labels = paste0("Training peptide\nlength: ", 
                                           c("6", "10 or less", "15 or less"))),
         et = factor(et, labels = c("Reduced alphabet", "Best performing reduced alphabet",
                                    "Reduced alphabet from literature", "Full alphabet"))) %>%
  mutate(len_range = factor(len_range, 
                            labels = paste0("Testing peptide\nlength: ", c("6 ", "7-10", "11-15", "16-25"))))

# Fig 1 all encodings sens/spec  ----------------------------------------

sesp_plot <- ggplot(amyloids_plot, aes(x = Sens_mean, y = Spec_mean, color = et, shape = et)) +
  geom_point() +
  scale_color_manual("", values = c("grey", "red", "blue", "green")) +
  scale_shape_manual("", values = c(1, 16, 15, 15)) +
  scale_x_continuous("Mean sensitivity") +
  scale_y_continuous("Mean specificity") +
  facet_grid(pos ~ len_range) +
  my_theme +
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
             aes(x = Sens_mean, y = Spec_mean, color = et))

#png("./pub_figures/sesp_plot.png", height = 648, width = 648)
#cairo_ps("./pub_figures/sesp_plot.eps", height = 9, width = 9)
# should be eps, but it's too big for overleaf
cairo_pdf("./pub_figures/sesp_plot.pdf", height = 9, width = 9)
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

cairo_ps("./pub_figures/AUC_boxplot.eps", height = 9, width = 9)
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

cairo_ps("./pub_figures/MCC_boxplot.eps", height = 9, width = 9)
#png("./pub_figures/MCC_boxplot.png", height = 648, width = 648)
MCC_boxplot
dev.off()