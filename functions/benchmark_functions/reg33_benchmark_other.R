# reg33 comes from http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054175
# it's named 33, but has 34 proteins
# add 1 do indices
# overlapping regions - what to do?

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)
require(hmeasure)
require(pbapply)

load("aa_groups.RData")
aa_groups <- string2list(aa_groups)


r33_raw <- readLines("./benchmark/reg33.txt") %>% matrix(ncol = 3, byrow = TRUE) %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  rename(name = X1, seq = X2, pos = X3)
# after analysis of table in PLoS article, r33_raw[c(5, 17), ] are wrongly annotated
# c("Serum Amyloid A", "Beta2-microglobulin")

hotspot_pos <- lapply(strsplit(r33_raw[, "pos"], " "), function(i) matrix(as.numeric(i), nrow = 2))
#check the lengths
r33_seqs <- strsplit(r33_raw[["seq"]], "")

r33_status <-lapply(1L:length(hotspot_pos), function(i) {
  # binary peptide (status of amyloidogenicity)
  bin_peptide <- rep(0, length(r33_seqs[[i]]))
  for (j in 1L:ncol(hotspot_pos[[i]])) {
    id_row <- hotspot_pos[[i]][, j]
    # important +1 to get proper position of 
    bin_peptide[id_row[1]:id_row[2] + 1] <- 1
  }
  bin_peptide
})

# low- and uppercase sequences to check with the source publication --------------------------------------------
r33_lu <-lapply(1L:length(hotspot_pos), function(i) {
  new_peptide <- r33_seqs[[i]]
  for (j in 1L:ncol(hotspot_pos[[i]])) {
    id_row <- hotspot_pos[[i]][, j]
    # important +1 to get proper position of 
    new_peptide[id_row[1]:id_row[2] + 1] <- tolower(new_peptide[id_row[1]:id_row[2] + 1])
  }
  new_peptide
})
write.fasta(r33_lu, r33_raw[["name"]], "./benchmark/reg33_lu.fasta")


# uppercase sequences for benchmark -----------------------------
write.fasta(r33_seqs, r33_raw[["name"]], "./benchmark/reg33.fasta")


# AmyLoad sequences ---------------------------------
raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]


# reg33 hot spots ------------------------------------

hs_lens <- lapply(hotspot_pos, function(i) i[2, ] - i[1, ])

reg33_summary <- data.frame(hs_len = unlist(hs_lens),
                            prot_len = rep(lengths(r33_seqs), lengths(hs_lens)),
                            prot_name = rep(r33_raw[["name"]], lengths(hs_lens)),
                            hs_n = factor(rep(lengths(hs_lens), lengths(hs_lens))))


# comparison of AmyLoad and reg33 ----------------------------------

# TRUE if the sequence from AmyloGram is also in reg33
reg33_AmyloGram <- sapply(seqs_list, function(i) {
  pasted_seq <- paste0(i, collapse = "")
  sum(grepl(pasted_seq, r33_raw[["seq"]])) > 0
})
# table(ets[reg33_AmyloGram])
# 467 sequences extracted from AmyLoad to train AmyloGram are in reg33 (355 nonamyloids and 112 amyloids)

status_reg33_AmyloGram <- sapply(seqs_list[reg33_AmyloGram], function(single_seq) {
  pasted_seq <- paste0(single_seq, collapse = "")
  seq_pos <- gsub(pasted_seq, paste0(rep("_", length(single_seq)), collapse = ""), r33_raw[["seq"]]) %>% 
    strsplit(split = "") %>% 
    lapply(function(i) i == "_")
  
  search_ids <- which(sapply(seq_pos, sum) > 1)
  lapply(search_ids, function(single_reg33_seq)
    r33_status[[single_reg33_seq]][seq_pos[[single_reg33_seq]]]
  )
})

decide_reg <- function(x) {
  ets <- unique(x)
  if(length(ets) == 1) {
    ets
  } else {
    2
  }
}

# duplicates: remove using data.frame(sort(table(reg33_al_comp[["seq"]]), decreasing = TRUE))

# hot spot fraction
calc_frac_hs <- function(x) {
  status <- as.numeric(unlist(strsplit(x, split = "")))
  sum(status)/length(status)
}

reg33_al_comp <- data.frame(seq_name = names(seqs_list[reg33_AmyloGram]),
                            seq = sapply(seqs_list[reg33_AmyloGram], paste0, collapse = ""),
                            reg33_status = sapply(status_reg33_AmyloGram, paste0, collapse = ""),
                            AmyLoad_et = ets[reg33_AmyloGram],
                            reg33_et = sapply(status_reg33_AmyloGram, decide_reg),
                            stringsAsFactors = FALSE
) %>% mutate(AmyLoad_et = factor(AmyLoad_et, labels = c("Non-amyloid", "Amyloid")),
             reg33_et = factor(reg33_et, labels = c("Non-amyloid", "Amyloid", "Mixed")),
             concordance = as.character(AmyLoad_et) == as.character(reg33_et)) %>% 
  slice(-c(16, 99, 272, 421)) %>% 
  group_by(seq) %>% 
  mutate(frac_hs = calc_frac_hs(reg33_status),
         AmyLoad_et_bin = ifelse(AmyLoad_et == "Amyloid", 1, 0),
         CI = abs(AmyLoad_et_bin - frac_hs)) %>% 
  ungroup %>% 
  select(-AmyLoad_et_bin, -frac_hs)


as.data.frame(table(reg33_al_comp[, c("AmyLoad_et", "reg33_et")]), responseName = "Count")

# training of AmyloGram ---------------------------

seq_hex <- lapply(r33_seqs, function(single_seq)
  seq2ngrams(single_seq, 6, a()[-1]) %>% 
    decode_ngrams() %>% 
    strsplit("") %>% 
    do.call(rbind, .)
) %>% 
  do.call(rbind, .) %>% 
  apply(2, tolower)

seq_hex_status <- lapply(r33_status, function(single_status)
  seq2ngrams(single_status, 6, C(0, 1)) %>% 
    decode_ngrams() %>% 
    strsplit("") %>% 
    do.call(rbind, .) %>% 
    as.numeric %>% 
    matrix(ncol = 6) %>% 
    rowMeans()
)

# reg33_preds <- make_classifier_whole_protein(tolower(t(sapply(seqs_list[!reg33_AmyloGram], 
#                                                               function(i) c(i, rep(NA, max(lengths(seqs_list[!reg33_AmyloGram])) - length(i)))))), 
#                                              ets[!reg33_AmyloGram], 
#                                              unname(lengths(seqs_list[!reg33_AmyloGram])), 
#                                              6, 
#                                              aa_groups[14592], 
#                                              seq_hex
# )
# 
# save(reg33_preds, file = "./results/reg33_preds.RData")
load("./results/reg33_preds.RData")

reg33_AmyloGram <- lapply(1L:length(seq_hex_status), function(i)
  data.frame(prot = i, status = seq_hex_status[[i]])) %>% 
  do.call(rbind, .) %>% 
  mutate(pred = reg33_preds) %>% 
  mutate(status_bin = floor(status))

reg33_AmyloGram_bench <- filter(reg33_AmyloGram, prot != 28)

#HMeasure(reg33_AmyloGram_bench[["status_bin"]], reg33_AmyloGram_bench[["pred"]], threshold = 0.5)[["metrics"]]

reg33_measures <- lapply(unique(reg33_AmyloGram_bench[["prot"]]), function(i) {
  dat <- reg33_AmyloGram_bench[reg33_AmyloGram_bench[["prot"]] == i, ]
  tryCatch(HMeasure(dat[["status_bin"]], dat[["pred"]])[["metrics"]], 
           warning = function(w) 
             list(HMeasure(dat[["status_bin"]], dat[["pred"]])[["metrics"]], w))
})

#check if scores were switched
reverted <- !sapply(reg33_measures, is.data.frame)

reg33_full_preds <- data.frame(prot = unique(reg33_AmyloGram_bench[["prot"]]),
                               reverted = reverted,
                               len = as.vector(table(reg33_AmyloGram_bench[["prot"]])) + 5,
                               hs_mean_len = sapply(hotspot_pos, function(i) mean(i[2, ] - i[1, ]))[-28],
                               hs_max_len = sapply(hotspot_pos, function(i) max(i[2, ] - i[1, ]))[-28],
                               do.call(rbind, lapply(1L:length(reverted), function(single_range_id) {
                                 if(reverted[single_range_id]) {
                                   reg33_measures[[single_range_id]][[1]]
                                 } else {
                                   reg33_measures[[single_range_id]]
                                 }
                               }))
) %>% mutate(len_d = cut(len, c(20, 51, 101, 170, 300, 600, 800)),
             hs_mean_len_d = cut(hs_mean_len, breaks = c(2, 6, 10, 15, 25, 50, 150)))


