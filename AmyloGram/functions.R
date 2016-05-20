library(signalHsmm)
library(biogram)
library(ranger)

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