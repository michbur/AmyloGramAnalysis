
ngram_status <- select(ngram_freq, decoded_name, diff_freq) %>% 
  mutate(status = ifelse(diff_freq > 0, "Amyloid", "Non-amyloid")) %>% 
  select(decoded_name, status) %>% 
  rename(ngram = decoded_name)

trans_df <- expand.grid(list(c(1L:6, "_"), c(1L:6, "_")), stringsAsFactors = TRUE) %>% 
  apply(1, paste0, collapse = "") %>% 
  sapply(function(single_transition) grepl(single_transition, levels(ngram_freq[["decoded_name"]]))) %>% 
  melt %>% 
  mutate(Var1 = factor(Var1, labels = levels(ngram_freq[["decoded_name"]]))) %>% 
  rename(ngram = Var1, trans = Var2) %>% 
  mutate(len = nchar(as.character(ngram))) %>% 
  filter(len > 1) %>%
  inner_join(ngram_status) %>% 
  group_by(trans, status) %>% 
  summarise(freq = mean(value)) %>% 
  mutate(start = substr(as.character(trans), 0, 1),
         end = substr(as.character(trans), 2, 2))

ggplot(trans_df, aes(x = status, y = freq, fill = freq, label = round(freq, 4))) +
  geom_bar(stat = "identity") +
  geom_text(color = "red") +
  facet_grid(start ~ end, labeller = label_both)
