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
  arrange(classifiers, len_range)

# benchmark table --------------------------------------------

read.csv("./results/benchmark_allpreds.csv")
