source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(cvTools)
require(ranger)
require(hmeasure)

ftraits <- choose_properties()

#time-consuming step commented out. Instead just load the results.
#aa_groups <- create_encodings(ftraits)
#save(aa_groups, file = "aa_groups.RData")
load("aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
               read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))

#for single aa_group around 85 sec
extracted_ngrams <- extract_ngrams(seqs_m, aa_groups)
#add to extracted_ngrams raw aa n-grams

#create folds for cv
ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]
seq_lengths <- unname(lengths(seqs_list))

all_folds <- create_all_folds(ets, seq_lengths)

cv_results <- do_cv(all_folds, extracted_ngrams, create_hv(seqs_m))

cv_summary <- summarize_cv_results(cv_results, get_reps_ids)

# negative control - full amino acid alphabet
# full_aa <- tolower(a()[-1]) %>% as.list
# names(full_aa) <- 1L:20
# full_aa <- list(full_aa)
# extracted_ngrams_full <- extract_ngrams(seqs_m, full_aa)
# cv_results_full <- do_cv(all_folds, extracted_ngrams_full, create_hv(seqs_m))
# save(cv_results_full, file = "./results/cv_results_full.RData")
load("./results/cv_results_full.RData")

