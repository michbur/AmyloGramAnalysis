source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)


load("aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_benchmark.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_benchmark.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 16
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_benchmark.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_benchmark.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

seq_lengths <- unname(lengths(seqs_list))

test_dat <- read.fasta("./benchmark/pep424_better_names.fasta")
test_dat_m <- tolower(t(sapply(test_dat, function(i)
  c(i, rep(NA, max(lengths(test_dat)) - length(i))))))

#aa_group must a list of length 1
make_classifier <- function(dat, ets, seq_lengths, max_len, aa_group, test_dat) {
  fdat <- extract_ngrams(dat[seq_lengths <= max_len, ], aa_group)[[1]]
  fets <- ets[seq_lengths <= max_len]
  
  test_bis <- test_features(fets, fdat, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  train_data <- data.frame(as.matrix(fdat[, imp_bigrams]), tar = factor(fets))
  
  model <- ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE)
  
  test_ngrams <- extract_ngrams(test_dat, aa_group)[[1]]
  
  test_lengths <- apply(test_dat, 1, function(i) sum(!is.na(i))) - 5
  
  preds <- cbind(predict(model, data.frame(as.matrix(test_ngrams)[, imp_bigrams]))[["predictions"]][, 2], 
                 unlist(lapply(1L:length(test_lengths), function(i) rep(i, test_lengths[i]))))
  
  preds %>% 
    data.frame %>% 
    rename(prob = X1, prot = X2) %>%
    group_by(prot) %>%
    # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
    # so we take maximum probabilities for all hexagrams belonging to the peptide
    summarise(prob = max(prob))
}

#time around 40 [s]
tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m)
system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))
