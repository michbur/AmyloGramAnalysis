source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")

compare_seq <- function(shorter, longer) {
  windows <- seq2ngrams(as.vector(longer), n = length(shorter), u = a()[-1]) %>% 
    decode_ngrams %>% 
    strsplit(split = "")
  
  sapply(windows, function(single_window) {
    unname(mean(shorter == single_window))
  }) %>% max
}

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)
require(hmeasure)
require(pbapply)

load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 16
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))))
ets <- ets[purified_seqs_id]

seq_lengths <- unname(lengths(seqs_list))

test_dat <- read.fasta("./benchmark/pep424_better_names.fasta", seqtype = "AA")
test_dat_m <- tolower(t(sapply(test_dat, function(i)
  c(i, rep(NA, max(lengths(test_dat)) - length(i))))))

pasted_seqs <- sapply(seqs_list, paste0, collapse = "")
pasted_test <- sapply(test_dat, paste0, collapse = "")

# not in the test set
nit <- sapply(pasted_seqs, function(i) sum(grepl(i, pasted_test, fixed = TRUE))) == 0

# time around 40 [s]
# tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m)
# system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))

# Get list of classifiers

full_aa <- tolower(a()[-1]) %>% as.list
names(full_aa) <- 1L:20
full_aa <- list(full_aa)

# window_pids <- pbsapply(test_dat, function(single_pep424) {
#   single_pep424 <- as.vector(single_pep424)
#   len_pep424 <- length(single_pep424)
#   
#   sapply(seqs_list, function(single_train) {
#     single_train <- as.vector(single_train)
#     len_train <- length(single_train)
#     
#     if(len_train < len_pep424) {
#       compare_seq(single_train, single_pep424)
#     } else {
#       compare_seq(single_pep424, single_train)
#     }
#   })
# })
# 
# save(window_pids, file = "./results/window_pids.RData")

load("./results/window_pids.RData")

class_list_test <- pblapply(c(0.7, 0.9, 1), function(threshold) {
  redundant_ids <- lapply(1L:ncol(window_pids), function(col_id) {
    dat <- window_pids[, col_id]
    which(dat > threshold)
  }) %>% unlist() %>% unique %>% sort
  
  if(length(redundant_ids) == 0) {
    pid_ids <- (1L:length(seqs_list))
  } else {
    pid_ids <- (1L:length(seqs_list))[-redundant_ids]
  }
  
  names(pid_ids) <- names(seqs_list[pid_ids])
  
  class14592 = make_classifier(seqs_m[names(pid_ids), ], ets[pid_ids],
                               seq_lengths[pid_ids], 6, aa_groups[14592], test_dat_m)
})

tmp <- do.call(cbind, class_list_test)
names(tmp) <- paste0("thr", c(0.7, 0.9, 1))

HMeasure(read.csv("results/benchmark_allpreds_raw.csv")[[1]], tmp)["metrics"]

pblapply(c(0.5, 0.7, 0.9, 1), function(threshold) {
  redundant_ids <- lapply(1L:ncol(window_pids), function(col_id) {
    dat <- window_pids[, col_id]
    which(dat > threshold)
  }) %>% unlist() %>% unique %>% sort
  
  if(length(redundant_ids) == 0) {
    pid_ids <- (1L:length(seqs_list))
  } else {
    pid_ids <- (1L:length(seqs_list))[-redundant_ids]
  }
  
  table(ets[pid_ids]) %>% 
    data.frame() %>% mutate(Freq = Freq/sum(Freq))
})

class_list <- pblapply(c(0.5, 0.7, 0.9, 1), function(threshold) {
  redundant_ids <- lapply(1L:ncol(window_pids), function(col_id) {
    dat <- window_pids[, col_id]
    which(dat > threshold)
  }) %>% unlist() %>% unique %>% sort
  
  if(length(redundant_ids) == 0) {
    pid_ids <- (1L:length(seqs_list))
  } else {
    pid_ids <- (1L:length(seqs_list))[-redundant_ids]
  }
  
  names(pid_ids) <- names(seqs_list[pid_ids])
  
  lapply(c(6, 10, 15), function(single_length)
    data.frame(class14592 = make_classifier(seqs_m[names(pid_ids), ], ets[pid_ids],
                                            seq_lengths[pid_ids], single_length, aa_groups[14592], test_dat_m),
               raw_aa = make_classifier(seqs_m[names(pid_ids), ], ets[pid_ids],
                                        seq_lengths[pid_ids], single_length, full_aa[1], test_dat_m)
    )
  )
})

save(class_list, file = "results/class_list_PID.RData")