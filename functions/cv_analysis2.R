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
require(e1071)

if(Sys.info()["nodename"] %in% c("phobos", "michal-XPS14"))
  pathway <- "/home/michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/amyloid_cv_results/"


extract_cv_results_single_file <- function(cv_results) {
  res <- do.call(rbind, lapply(1L:length(cv_results), function(single_enc_id) {
    all_reps_res <- do.call(rbind, lapply(cv_results[[single_enc_id]], function(single_rep) {
      do.call(rbind, lapply(single_rep, function(single_fold_res) {
        single_fold_res[[1]]
      })) %>% mutate(MCC = calc_mcc(TP, TN, FP, FN),
                     fold = unlist(lapply(1L:5, rep, 4)),
                     enc = single_enc_id) %>%
        select(enc, fold, len_range, AUC, MCC, Sens, Spec)
    }))
  }))
  res
}


all_summaries <- lapply(1L:20, function(single_partition_id) {
  single_part_res <- do.call(rbind, pblapply(1L:45, function(single_replicate_id) {
    load(paste0(pathway, "results/cv_results_full_", single_replicate_id, "_", 
                single_partition_id,".Rdata"))
    data.frame(replicate = floor((single_replicate_id - 1)/5) + 1, 
               extract_cv_results_single_file(cv_results))
  }))
  group_by(single_part_res, replicate, len_range, enc) %>%
    summarize_each(funs(mean, sd), AUC, MCC, Sens, Spec) %>%
    mutate(partition  = single_partition_id)
})
