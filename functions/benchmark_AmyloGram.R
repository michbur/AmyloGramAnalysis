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

class_list <- lapply(c(6, 10, 15), function(single_length)
  data.frame(class15608 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[15608], test_dat_m),
             class3276 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[3276], test_dat_m),
             class16792 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[16792], test_dat_m),
             class12626 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[12626], test_dat_m),
             class14592 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14592], test_dat_m),
             class14596 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14596], test_dat_m),
             class14777 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14777], test_dat_m),
             class10880 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[10880], test_dat_m))
)
#time around 40 [s]
# tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m)
# system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))

save(class_list, file = "./results/class_list.RData")

# dat <- cbind(read.csv("./results/benchmark_otherpreds.csv"),
#              class15608,
#              class3276,
#              class16792,
#              class12626,
#              class14592,
#              class14596,
#              class14777,
#              class10880)
# 
# HMeasure(dat[[1]], dat[-1])[["metrics"]] %>%
#   mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
#   select(classifier, AUC, MCC, Sens, Spec)
