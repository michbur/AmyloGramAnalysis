#' Explain replicates of cross-validation
#'
#' Since the integration with WCSS, the exact structure of cross-validation results 
#' might be obstruse, so this function additional information about the cross-validation 
#' replicates.
#' 
#' @return List of maximum lengths of sequences in training set.

get_reps_ids <- function() {
  reps_ids <- lapply(c(6, 10, 15), function(constant_pos) {
    lapply(c(6, 10, 15), function(constant_neg) {
      lapply(1L:5, function(dummy) {
        c(pos = constant_pos, neg = constant_neg)
      })
    })
  }) %>% unlist %>% 
  matrix(ncol = 2, byrow = TRUE) %>%
  data.frame %>%
  rename(pos = X1, neg = X2)

  reps_ids[sapply(1L:nrow(reps_ids), rep, 4), ]
}

#helper function calculating MCC
calc_mcc <- function(TP, TN, FP, FN)
  (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))


#' Summarize results of cross-validation
#'
#' Calculates performance measures of cross-valiation results.
#' @param cv_results as returned by do_cv
#' @param get_reps_ids function returning the lengths of negative and positive
#' sequences in the training set. For example get_reps_ids.
#' 
#' @return Data frame of performance measures - their mean values and standard 
#' deviations.
summarize_cv_results <- function(cv_results, get_reps_ids) {
  res <- do.call(rbind, lapply(1L:length(cv_results), function(single_enc_id) {
    all_reps_res <- lapply(cv_results[[single_enc_id]], function(single_rep) {
      do.call(rbind, lapply(single_rep, function(single_fold_res) {
        single_fold_res[[1]]
      })) %>% mutate(MCC = calc_mcc(TP, TN, FP, FN)) %>%
        select(len_range, AUC, MCC, Sens, Spec)
    }) %>% do.call(rbind, .) %>% 
      cbind(., get_reps_ids()) %>% 
      group_by(len_range, pos, neg) %>%
      summarize_each(funs(mean, sd)) %>%
      cbind(., enc = single_enc_id)
  }))

  # dplyr changes order of len_ranger, need to get it back
  res[["len_range"]] <- factor(res[["len_range"]], levels(res[["len_range"]])[c(4, 3, 1, 2)])
  res[["pos"]] <- as.factor(res[["pos"]])
  res[["neg"]] <- as.factor(res[["neg"]])
  res
}



# levels(cv_summary[["len_range"]]) <- paste0("Test peptide length: ", levels(cv_summary[["len_range"]]))
# library(ggplot2)
# ggplot(cv_summary, aes(x = pos, y = neg, fill = mauc, label = round(mauc, 4))) +
#   geom_tile() +
#   geom_text(color = "red") +
#   facet_wrap(~ len_range) +
#   scale_x_discrete("Positive training set") +
#   scale_y_discrete("Positive negative set")


