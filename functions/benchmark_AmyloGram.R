source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")


require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)
require(hmeasure)
require(pbapply)

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


# time around 40 [s]
# tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m)
# system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))

# Get list of classifiers

full_aa <- tolower(a()[-1]) %>% as.list
names(full_aa) <- 1L:20
full_aa <- list(full_aa)

class_list <- list(
  class14592 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14592], test_dat_m),
  class14596 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14596], test_dat_m),
  class14533 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14533], test_dat_m),
  class18297 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[18297], test_dat_m),
  class6704 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[6704], test_dat_m),
  raw_aa = make_classifier(seqs_m, ets, seq_lengths, single_length, full_aa[1], test_dat_m)
)

save(class_list, file = "./results/class_list_best5_full.RData")

# learn_lengths <- c(6, 10, 15)
# load("./results/class_list_best5.RData")
# 
# 
# dat <- cbind(read.csv("./results/benchmark_otherpreds.csv"),
#              do.call(cbind, lapply(1L:3, function(i) {
#                single_preds_df <- class_list[[i]]
#                colnames(single_preds_df) <- paste0(colnames(single_preds_df), "_", learn_lengths[i])
#                single_preds_df
#              })))
# 
# HMeasure(dat[[1]], dat[-1])[["metrics"]] %>%
#   mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
#   select(classifier, AUC, MCC, Sens, Spec) 
