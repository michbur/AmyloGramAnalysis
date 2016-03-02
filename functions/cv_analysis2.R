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

if(Sys.info()["nodename"] == "tobit" )
  pathway <- "~/Dropbox/amyloid_cv_results/"

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


get_reps_ids <- function() 
  lapply(c(6, 10, 15), function(constant_pos) {
    lapply(c(6, 10, 15), function(constant_neg) {
      c(pos = constant_pos, neg = constant_neg)
    })
  }) %>% unlist %>% 
  matrix(ncol = 2, byrow = TRUE) %>%
  data.frame %>%
  rename(pos = X1, neg = X2)

all_summaries <- do.call(rbind, lapply(1L:20, function(single_partition_id) {
  single_part_res <- do.call(rbind, pblapply(1L:45, function(single_replicate_id) {
    load(paste0(pathway, "results/cv_results_full_", single_replicate_id, "_", 
                single_partition_id,".Rdata"))
    data.frame(replicate = floor((single_replicate_id - 1)/5) + 1, 
               extract_cv_results_single_file(cv_results))
  }))
  group_by(single_part_res, replicate, len_range, enc) %>%
    summarize_each(funs(mean, sd), AUC, MCC, Sens, Spec) %>%
    mutate(partition  = single_partition_id)
}))

enc_part_id_adj <- c(0, cumsum(rep(927, 19)))

all_summaries[["len_range"]] <- factor(all_summaries[["len_range"]], levels(all_summaries[["len_range"]])[c(4, 3, 1, 2)])
all_summaries <- mutate(all_summaries, 
                        pos = get_reps_ids()[replicate, "pos"],
                        neg = get_reps_ids()[replicate, "neg"],
                        enc_adj = enc + enc_part_id_adj[partition]) %>%
  select(len_range, enc_adj, pos, neg, AUC_mean, MCC_mean, Sens_mean, 
         Spec_mean, AUC_sd, MCC_sd, Sens_sd, Spec_sd, enc, partition)

write.csv(all_summaries, file = "./results/all_summaries.csv", row.names = FALSE)

best_sens <- group_by(all_summaries, pos, neg, len_range) %>%
  filter(Sens_mean == max(Sens_mean, na.rm = TRUE)) %>%
  arrange(desc(Sens_mean)) %>%
  select(pos, neg, len_range, AUC_mean, Sens_mean, Spec_mean, enc_adj) %>%
  filter(AUC_mean == max(AUC_mean)) %>%
  mutate(type = "Best sensitivity")

best_spec <- group_by(all_summaries, pos, neg, len_range) %>%
  filter(Spec_mean == max(Spec_mean, na.rm = TRUE)) %>%
  arrange(desc(Spec_mean)) %>%
  select(pos, neg, len_range, AUC_mean, Sens_mean, Spec_mean, enc_adj) %>%
  filter(AUC_mean == max(AUC_mean)) %>%
  mutate(type = "Best specificity")

bests <- data.frame(rbind(best_sens, best_spec))

group_by(bests, type, enc_adj) %>%
  summarize(occ = length(enc_adj)) %>%
  filter(occ > 1) %>%
  data.frame()

library(ggplot2)
ggplot(bests, aes(x = factor(pos), y = factor(neg), fill = Sens_mean, label = round(Sens_mean, 2))) +
  geom_tile() +
  geom_text(color = "red", size = 5) +
  facet_grid(type ~ len_range)

png("auc_plot.png", height = 1024, width = 1024)
ggplot(bests, aes(x = factor(pos), y = factor(neg), fill = AUC_mean, 
                  label = paste0("Se: ", formatC(Sens_mean, digits = 2, format = "f"), 
                                 "\n",
                                 "Sp: ", formatC(Spec_mean, digits = 2, format = "f")))) +
  geom_tile(color = "orange") +
  geom_text(color = "red", size = 5) +
  facet_grid(type ~ len_range)
dev.off()