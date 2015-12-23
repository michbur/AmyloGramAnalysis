reps_ids <- lapply(c(6, 10, 15), function(constant_pos) {
  lapply(c(6, 10, 15), function(constant_neg) {
    lapply(1L:10, function(dummy) {
      c(pos = constant_pos, neg = constant_neg)
    })
  })
}) %>% unlist %>% 
  matrix(ncol = 2, byrow = TRUE) %>%
  data.frame %>%
  rename(pos = X1, neg = X2)



single_enc <- cv_results_full[[1]]

all_reps_res <- lapply(single_enc, function(single_rep) 
  lapply(single_rep, function(single_fold_res) {
    aggr_preds <- single_fold_res[[1]]  %>% 
      data.frame %>% 
      rename(prob = X1, tar = X2, prot = X3) %>% 
      mutate(prot = paste0(tar, "_", prot)) %>%
      group_by(prot) %>%
      # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
      # so we take maximum probabilities for all hexagrams belonging to the peptide
      summarise(prob = max(prob), tar = unique(tar), len = 5 + length(prot)) %>%
      mutate(len_range = cut(len, include.lowest = TRUE, breaks = c(5, 6, 10, 15, max(len))))

    data.frame(len_range = levels(aggr_preds[["len_range"]]), 
               n = as.vector(table(aggr_preds[["len_range"]])),
               do.call(rbind, lapply(levels(aggr_preds[["len_range"]]), function(single_range) {
                 dat <- aggr_preds[aggr_preds[["len_range"]] == single_range, ]
                 HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]]
               })))
  })
)


data.frame(reps_ids, all_reps_res) %>% 
  select(pos, neg, AUC, Sens, Spec) %>%
  group_by(pos, neg) %>%
  summarise(mAUC = mean(AUC), mSens = mean(Sens), mSpec= mean(Spec))
