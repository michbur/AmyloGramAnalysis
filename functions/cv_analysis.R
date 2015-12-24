#' Explain replicates of cross-validation
#'
#' Since the integration with WCSS, the exact structure of cross-validation results 
#' might be obstruse, so this function additional information about the cross-validation 
#' replicates.
#' 
#' @return List of maximum lengths of sequences in training set.

get_reps_ids <- function()
  lapply(c(6, 10, 15), function(constant_pos) {
    lapply(c(6, 10, 15), function(constant_neg) {
      lapply(1L:10, function(dummy) {
        c(pos = constant_pos, neg = constant_neg)
      })
    })
  }) %>% unlist %>% 
  matrix(ncol = 2, byrow = TRUE) %>%
  data.frame %>%
  rename(pos = X1, neg = X2)

reps_ids <- get_reps_ids()


summarize_cv_results <- function(cv_results, reps_ids) {
  res <- do.call(rbind, lapply(1L:length(cv_results), function(single_enc_id) {
    all_reps_res <- lapply(cv_results[[single_enc_id]], function(single_rep) {
      all_folds <- do.call(rbind, lapply(single_rep, function(single_fold_res) {
        aggr_preds <- single_fold_res[[1]]  %>% 
          data.frame %>% 
          rename(prob = X1, tar = X2, prot = X3) %>% 
          mutate(prot = paste0(tar, "_", prot)) %>%
          group_by(prot) %>%
          # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
          # so we take maximum probabilities for all hexagrams belonging to the peptide
          summarise(prob = max(prob), tar = unique(tar), len = 5 + length(prot)) %>%
          mutate(len_range = cut(len, include.lowest = TRUE, breaks = c(5, 6, 10, 15, max(len))))
        
        # len_range
        levels(aggr_preds[["len_range"]])[4] <- "(15,max]"
        
        data.frame(len_range = levels(aggr_preds[["len_range"]]), 
                   n = as.vector(table(aggr_preds[["len_range"]])),
                   do.call(rbind, lapply(levels(aggr_preds[["len_range"]]), function(single_range) {
                     dat <- aggr_preds[aggr_preds[["len_range"]] == single_range, ]
                     HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]]
                   })))
      }))
      
      group_by(all_folds, len_range) %>%
        summarize(mn = mean(n), mauc = mean(AUC), mspec = mean(Spec), msens = mean(Sens))
    })
    
    do.call(rbind, suppressWarnings(lapply(1L:length(all_reps_res), function(id) {
      data.frame(all_reps_res[[id]], reps_ids[id, ])
    }))) %>% 
      group_by(len_range, pos, neg) %>%
      summarize(mn = mean(mn), mauc = mean(mauc), mspec = mean(mspec), msens = mean(msens)) %>%
      cbind(., enc = single_enc_id)
  }))
  
  # dplyr changes order of len_ranger, need to get it back
  res[["len_range"]] <- factor(res[["len_range"]], levels(res[["len_range"]])[c(3, 4, 1, 2)])
  res[["pos"]] <- as.factor(res[["pos"]])
  res[["neg"]] <- as.factor(res[["neg"]])
  res
}

cv_summary <- summarize_cv_results(cv_results_full, reps_ids)

levels(cv_summary[["len_range"]]) <- paste0("Test peptide length: ", levels(cv_summary[["len_range"]]))
library(ggplot2)
ggplot(cv_summary, aes(x = pos, y = neg, fill = mauc, label = round(mauc, 4))) +
  geom_tile() +
  geom_text(color = "red") +
  facet_wrap(~ len_range) +
  scale_x_discrete("Positive training set") +
  scale_y_discrete("Positive negative set")

