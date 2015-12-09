#' Perform cross-validation
#'
#' Performs cross-validation using the n-grams extracted from 
#' amylyoids and non-amyloids.
#' @param extracted_ngrams a list of n-grams as fiven by \code{\link{extract_ngrams}}.
#' @param ets a vector of etiquettes (0 or 1).
#' @param ets a vector of lengths of sequences.
#' 
#' @note The classifier of choice is the \code{\link{[ranger]ranger}} implementation of 
#' random forests algorithm.
#' 
#' @return List of length equal to the number of repetitions of cross-validation.
#' Each element of the list is a list with the number of elements equal to the number 
#' of encodings. Each element of this list contains a list with two elements:
#' a) data frame with the prediction, real value of target and id of the protein,
#' b) list of important n-grams as found by \code{\link[biogram]{test_features}}.

do_cv <- function(extracted_ngrams, ets, seq_lengths) {
  
  pos_train <- which(ets == 1 & seq_lengths < 7)
  neg_train <- which(ets == 0 & seq_lengths < 16)
  
  pos_test <- which(ets == 1 & seq_lengths >= 7)
  neg_test <- which(ets == 0 & seq_lengths >= 16)
  
  lapply(1L:2, function(dummy) {
    
    fold_list <- lapply(list(pos_train, pos_test, neg_train, neg_test), function(single_n) {
      folded <- cvFolds(length(single_n), K = 5)
      data.frame(id = single_n[folded[["subsets"]]], which = folded[["which"]])
    })
    
    lapply(1L:5, function(fold) {
      lapply(extracted_ngrams, function(encoded_group) {
        
        #training data
        train_pos <- do.call(rbind, encoded_group[fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"]])
        train_neg <- do.call(rbind, encoded_group[fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"]])
        
        test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                                  rbind(train_pos, train_neg), adjust = NULL)
        imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
        
        train_data <- data.frame(rbind(train_pos, train_neg)[, imp_bigrams],
                                 tar = as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
        model_cv <- ranger(tar ~ ., train_data, write.forest = TRUE)
        
        #test data
        test_pos <- rbind(do.call(rbind, encoded_group[fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]]),
                          do.call(rbind, encoded_group[fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]]))
        test_neg <- rbind(do.call(rbind, encoded_group[fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]]), 
                          do.call(rbind, encoded_group[fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]]))
        
        #number of n-grams from protein
        ngram_prots_pos <- c(sapply(encoded_group[fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]], nrow),
                             sapply(encoded_group[fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]], nrow))
        ngram_prots_neg <- c(sapply(encoded_group[fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]], nrow),
                             sapply(encoded_group[fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]], nrow))
        
        preds <- cbind(as.numeric(as.character(predict(model_cv, data.frame(rbind(test_pos, test_neg)[, imp_bigrams]))[["predictions"]])), 
                       c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))),
                       c(unlist(lapply(1L:length(ngram_prots_pos), function(prot_id)
                         rep(prot_id, ngram_prots_pos[prot_id]))), 
                         unlist(lapply(1L:length(ngram_prots_neg), function(prot_id)
                           rep(prot_id, ngram_prots_neg[prot_id])))))
        
        # is randomForest consistent with ranger?
        # model_rf <- randomForest(tar ~ ., train_data)
        # predict(model_rf, data.frame(rbind(test_pos, test_neg)[, imp_bigrams])) == 
        #   predict(model_cv, data.frame(rbind(test_pos, test_neg)[, imp_bigrams]))[["predictions"]]
        
        list(preds, imp_bigrams)
      })
    })
  })
}