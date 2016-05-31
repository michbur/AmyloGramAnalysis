library(seqinr)
library(dplyr)
library(biogram)

preds <- read.csv("./results/AmyloGram_pep424preds.csv")
pep424 <- read.fasta("./benchmark/pep424_better_names.fasta", seqtype = "AA")

pep424_m <- tolower(t(sapply(pep424, function(i)
  c(i, rep(NA, max(lengths(pep424)) - length(i)))
)))

ngram_freq <- read.csv("./results/ngram_freq.csv")

imp_ngrams <- ngram_freq[["decoded_name"]][ngram_freq[["association"]] != "Not found"] %>% as.character %>% code_ngrams

pep424_ngrams <- as.matrix(count_specified(degenerate(pep424_m, aa_groups[[best_enc]]), imp_ngrams))
colnames(pep424_ngrams) <- decode_ngrams(colnames(pep424_ngrams))
dat <- cbind(preds, pep424_ngrams) %>% 
  mutate(Name = as.character(Name),
         len = lengths(pep424)) %>% 
  arrange(Probability) %>% 
  mutate(Name = factor(Name, levels = Name))

library(ggplot2)

p <- ggplot(dat, aes(x = Name, y = Probability, color = status_bin)) +
  geom_point() +
  geom_hline(yintercept = 0.5) +
  my_theme +
  theme(panel.grid.major.x = element_line(color = NA))

library(plotly)
ggplotly(p, tooltip = c("x", "y"))

mdat <- select(dat, -Name) %>% 
  melt(id.vars = c("Probability", "status_bin", "len"), variable.name = "ngram") %>% 
  mutate(n = nchar(as.character(ngram))) %>% 
  group_by(ngram, status_bin) %>% 
  mutate(freq = value/(len - n + 1),
         cvalue = cumsum(freq),
         amyloid_ngram = ifelse(ngram %in% c("2", "2___2", "22"), "Non-amyloidogenic ngram", "Amyloidogenic ngram"),
         amyloid_status = ifelse(status_bin, "Amyloidogenic peptide", "Non-amyloidogenic peptide")
  )

p <- ggplot(filter(mdat, !(ngram %in% c("2", "3"))), 
            aes(x = cvalue, y = Probability, color = ngram)) +
  geom_step() +
  #geom_point() +
  facet_grid(amyloid_ngram ~ amyloid_status) +
  scale_x_continuous("Cumulative frequency") +
  my_theme

ggplotly(p, tooltip = c("x", "y"))

