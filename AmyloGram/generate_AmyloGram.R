library(seqinr)
library(biogram)
library(ranger)

source("./functions/encode_amyloids.R")
source("./functions/make_classifier.R")

load("aa_groups.RData")
best_encoding <- string2list(aa_groups)[[14592]]

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_benchmark.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_benchmark.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 16
seqs_list <- raw_seqs_list[purified_seqs_id]

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_benchmark.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_benchmark.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

test_dat <- read.fasta("./benchmark/pep424_better_names.fasta")
test_dat_m <- tolower(t(sapply(test_dat, function(i)
  c(i, rep(NA, max(lengths(test_dat)) - length(i))))))

preds <- make_classifier(seqs_m, ets, seq_lengths, 10, aa_groups[14592], test_dat_m)

AmyloGram_res <- data.frame(pred = preds,
                            status_bin = read.csv("./results/benchmark_otherpreds.csv")[[1]])


make_AmyloGram <- function(seqs_list, ets, max_len, aa_group) {
  seq_lengths <- lengths(seqs_list)

  seqs_m <- tolower(t(sapply(seqs_list[seq_lengths <= max_len], function(i)
    c(i, rep(NA, max_len - length(i)))
    )))
  
  gl <- do.call(rbind, lapply(1L:nrow(seqs_m), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seqs_m[i, ][!is.na(seqs_m[i, ])], 6, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  }))

  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                                    ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                                    seq = degenerate(gl[, -7], aa_group),
                                                    u = as.character(1L:length(aa_group))))
  fdat <- bitrigrams > 0
  storage.mode(fdat) <- "integer"
  fets_raw <- ets[seq_lengths <= max_len]
  flens <- seq_lengths[seq_lengths <= max_len] - 5
  fets <- unlist(lapply(1L:length(flens), function(i) rep(fets_raw[i], flens[i])))
  
  test_bis <- test_features(fets, fdat, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  train_data <- data.frame(as.matrix(fdat[, imp_bigrams]), tar = factor(fets))
  
  list(rf = ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE),
       imp_features = imp_bigrams,
       enc = aa_group)
}

AmyloGram_model <- make_AmyloGram(seqs_list, ets, 10, best_encoding)


predict_AmyloGram <- function(model, seqs_list) {
  seqs_m <- tolower(t(sapply(seqs_list, function(i)
    c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))
  
  gl <- do.call(rbind, lapply(1L:nrow(seqs_m), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seqs_m[i, ][!is.na(seqs_m[i, ])], 6, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  }))
  
  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           seq = degenerate(gl[, -7], model[["enc"]]),
                                           u = as.character(1L:length(model[["enc"]]))))
  
  test_ngrams <- bitrigrams > 0
  storage.mode(test_ngrams) <- "integer"
  
  test_lengths <- lengths(seqs_list) - 5
  
  preds <- data.frame(prob = predict(model[["rf"]], data.frame(test_ngrams[, model[["imp_features"]]]))[["predictions"]][, 2], 
                      prot = unlist(lapply(1L:length(test_lengths), function(i) rep(i, test_lengths[i])))
  )
  
  data.frame(seq_name = names(seqs_list),
             prob = vapply(unique(preds[["prot"]]), function(single_prot)
    max(preds[preds[["prot"]] == single_prot, "prob"]),
    0)
  )
}

rep424_preds <- predict_AmyloGram(AmyloGram_model, test_dat)
rep424_preds[["status_bin"]] <- as.numeric(read.csv("./results/benchmark_otherpreds.csv")[[1]])

balance_spec_sens <- function(x, cutoffs) {
  do.call(rbind, lapply(cutoffs, function(cutoff) {
    x[["pred_bin"]] <- as.numeric(x[["prob"]] > cutoff)
    conf_mat <- as.data.frame(table(pred = x[["pred_bin"]], et = x[["status_bin"]]), responseName = "count")
    conf_mat[["count"]] <- as.numeric(conf_mat[["count"]])
    
    data.frame(Cutoff = cutoff,
               Sensitivity = conf_mat[4, "count"]/(conf_mat[4, "count"] + conf_mat[3, "count"]), 
               Specificity = conf_mat[1, "count"]/(conf_mat[1, "count"] + conf_mat[2, "count"]),
               MCC = calc_mcc(conf_mat[4, "count"], conf_mat[1, "count"], conf_mat[2, "count"], conf_mat[3, "count"]))
  }))
}

spec_sens <- balance_spec_sens(rep424_preds, 1L:99/100)

make_decision(predict_AmyloGram(AmyloGram_model, test_dat[1L:10]), 0.5)

make_decision <- function(x, cutoff) 
  data.frame(x, amyloidogenic = x[["prob"]] > cutoff)

save(spec_sens, AmyloGram_model, file = "./AmyloGram/AmyloGram.RData")