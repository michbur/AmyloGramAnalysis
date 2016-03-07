format_table <- function(x, caption, label, range) {
  bold_max <- function(x) {
    nx <- as.numeric(x)
    x[nx == max(nx)] <- paste0("\\textbf{", x[nx == max(nx)], "}")
    x
  }
  
  pretty_format <- function(x)
    formatC(x, digits = 4, format = "f")
  
  tab_bold <- x[range] %>%
    mutate_each(funs(pretty_format)) %>% 
    mutate_each(funs(bold_max))
  
  rws <- seq(1, nrow(tab_bold) - 1, by = 2)
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  
  tab_full <- cbind(x[-range], tab_bold)
  
  res <- print(xtable(tab_full, caption = caption, label = label), include.rownames = FALSE, booktabs = TRUE,
               add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE, 
               caption.placement = "top", sanitize.text.function = identity, 
               sanitize.rownames.function = identity)
  res
}

# cross-validation summary table ----------------------------

pub_full_summary_table <- rbind(select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>%
                                  group_by(pos, len_range) %>%
                                  summarise_each(funs(mean)) %>%
                                  mutate(classifiers = "All reduced alphabets"),
                                select(amyloids, pos, len_range, AUC_mean, MCC_mean, Sens_mean, Spec_mean, enc_adj) %>%
                                  filter(enc_adj == best_enc) %>%
                                  select(-enc_adj) %>%
                                  mutate(classifiers = "Best reduced alphabet"),
                                select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
                                  filter(enc_adj %in% c(1L:2)) %>%
                                  select(-enc_adj) %>%
                                  group_by(pos, len_range) %>%
                                  summarise_each(funs(mean)) %>%
                                  mutate(classifiers = "Reduced alphabet from literature"),
                                select(full_alphabet, pos, len_range, AUC_mean, MCC_mean, Sens_mean, Spec_mean) %>%
                                  mutate(classifiers = "Full alphabet")
)

colnames(pub_full_summary_table) <- c("pos", "len_range", "AUC_mean", "MCC_mean", "Sens_mean", "Spec_mean", 
                                      "classifiers")



mutate(pub_full_summary_table, 
       len_range = factor(len_range, levels = levels(len_range)[c(3L:4, 1L:2)])) %>%
  group_by(len_range, classifiers) %>%
  summarise_each(funs(mean)) %>%
  ungroup %>%
  arrange(classifiers, len_range) %>%
  select(-pos) %>%
  rename("Length of peptides" = len_range, Classifier = classifiers, "Mean AUC" = AUC_mean,
         "Mean MCC" = MCC_mean, "Mean sensitivity" = Sens_mean, "Mean specificity" = Spec_mean) %>%
  select(2, 1, 3L:6)
  format_table("Results of cross-validation.", "table:cv_summary", 3L:6) %>%
  cat



# benchmark table --------------------------------------------


