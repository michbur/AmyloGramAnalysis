all_prots

pep1 <- c("a", "a", "a", "a", "a", "a", "a")
pep2 <- c("a", "a", "a", "a", "a", "a")
calc_hamming(pep1, pep2) 

calc_hamming <- function(x, y) {
  if(length(x) < length(y)) {
    x <- c(x, rep("-", length(y) - length(x)))
  } else {
    y <- c(y, rep("-", length(x) - length(y)))
  }
  
  hamming.distance(x, y) 
}

hd1 <- pblapply(test_dat, function(single_pep424) {
  single_pep424 <- as.vector(single_pep424)

  which(sapply(seqs_list, function(single_train) {
    single_train <- as.vector(single_train)

    calc_hamming(single_pep424, single_train) == 1
  }))
})


et_pep424 <- read.csv("results/benchmark_allpreds_raw.csv")[[1]]

hd1_prots <- lapply(1L:length(hd1), function(i) {
  train_ids <- hd1[[i]]
    
  list(pep424_prot = as.vector(test_dat[[i]]),
       pep424_label = et_pep424[i],
       #train_prots = seqs_list[train_ids],
       train_et = data.frame(table(factor(as.logical(ets[train_ids]), levels = c(FALSE, TRUE)))))
})

library(ggplot2)

h1_dat <- lapply(hd1_prots[lengths(hd1) != 0], function(i) {
  data.frame(prot = paste0(i[["pep424_prot"]], collapse = ""),
             et = i[["pep424_label"]],
             len = length(i[["pep424_prot"]]),
             i[["train_et"]])
}) %>% 
  do.call(rbind, .) %>% 
  group_by(prot) %>% 
  mutate(Freq = Freq/sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Var1 = factor(Var1, labels = c("Non-amyloid", "amyloid")),
         et = factor(et, labels = c("Non-amyloid peptides (pep424 data set)", 
                                    "Amyloid peptides (pep424 data set)")))

prot_order <- filter(h1_dat, Var1 == "Non-amyloid") %>% 
  arrange(desc(Freq)) %>% 
  select(prot) %>%
  unlist() %>% 
  as.character()

cairo_ps("hamming1_pep424.eps", height = 8.1, width = 5.3)
ggplot(mutate(h1_dat, prot = factor(prot, levels = prot_order)), 
       aes(x = prot, y = Freq, fill = Var1)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~ et, ncol = 1, scales = "free_y") +
  coord_flip() +
  scale_y_continuous("Fraction of peptides in Hamming distance 1") +
  scale_x_discrete("") +
  scale_fill_manual("Status of peptide in the training data set", values = c("blue", "red")) +
  my_theme
dev.off()


hd1_train <- pblapply(seqs_list, function(single_pep424) {
  single_pep424 <- as.vector(single_pep424)
  
  which(sapply(seqs_list, function(single_train) {
    single_train <- as.vector(single_train)
    
    calc_hamming(single_pep424, single_train) == 1
  }))
})