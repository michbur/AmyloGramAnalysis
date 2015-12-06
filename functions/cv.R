#ranger

fold_res <- lapply(c(6, 10, 15), function(constant) {
  
  pos_train <- which(ets == 1 & seq_lengths < 7)
  neg_train <- which(ets == 0 & seq_lengths < 16)
  
  
  pos_constant <- apply(pos_data, 1, function(i) length(i) - sum(is.na(i))) <= constant
  neg_constant <- apply(neg_data, 1, function(i) length(i) - sum(is.na(i))) > constant
  
  lapply(1L:8, function(dummy) {
    fold_list <- lapply(list(pos_constant, !pos_constant, neg_constant, !neg_constant), function(single_n) {
      folded <- cvFolds(sum(single_n), K = 5)
      cbind(id = which(single_n)[folded[["subsets"]]], which = folded[["which"]])
    })
    
    lapply(1L:5, function(fold) {
      lapply(1L:length(aa_groups), function(group_id) {
        train_pos <- do.call(rbind, coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"]])
        train_neg <- do.call(rbind, coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"]])
        
        test_pos <- rbind(do.call(rbind, coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]]),
                          do.call(rbind, coded_pos[[group_id]][fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]]))
        test_neg <- rbind(do.call(rbind, coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]]), 
                          do.call(rbind, coded_neg[[group_id]][fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]]))
        
        test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                                  rbind(train_pos, train_neg), adjust = NULL)
        imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
        
        
        model_cv <- ranger(x = rbind(train_pos, train_neg)[, imp_bigrams], 
                                 as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
        
        #number of n-grams from protein
        ngram_prots_pos <- c(sapply(coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]], nrow),
                             sapply(coded_pos[[group_id]][fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]], nrow))
        ngram_prots_neg <- c(sapply(coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]], nrow),
                             sapply(coded_neg[[group_id]][fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]], nrow))
        
        preds <- cbind(predict(model_cv, rbind(test_pos, test_neg), type = "prob")[, 2], 
                       c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))),
                       c(unlist(lapply(1L:length(ngram_prots_pos), function(prot_id)
                         rep(prot_id, ngram_prots_pos[prot_id]))), 
                         unlist(lapply(1L:length(ngram_prots_neg), function(prot_id)
                           rep(prot_id, ngram_prots_neg[prot_id])))))
        
        
        list(preds, imp_bigrams)
      })
    })
  })
})