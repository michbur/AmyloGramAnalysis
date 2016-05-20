library(seqinr)
library(biogram)
library(ranger)

source("./functions/encode_amyloids.R")
source("./functions/make_classifier.R")
source("./AmyloGram/functions.R")

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

AmyloGram_model <- make_AmyloGram(seqs_list, ets, 10, best_encoding)

rep424_preds <- predict_AmyloGram(AmyloGram_model, test_dat)
rep424_preds[["status_bin"]] <- 

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
