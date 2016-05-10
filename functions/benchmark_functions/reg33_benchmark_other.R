# reg33 comes from http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054175
# it's named 33, but has 34 proteins
# add 1 do indices
# overlapping regions - what to do?
#read reg33 file

library(dplyr)
library(seqinr)


r33_raw <- readLines("./benchmark/reg33.txt") %>% matrix(ncol = 3, byrow = TRUE) %>% 
  data.frame(stringsAsFactors = FALSE) %>% 
  rename(name = X1, seq = X2, pos = X3)
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

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

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

reg33_al_comp <- data.frame(seq_name = names(seqs_list[reg33_AmyloGram]),
                            seq = sapply(seqs_list[reg33_AmyloGram], paste0, collapse = ""),
                            reg33_status = sapply(status_reg33_AmyloGram, paste0, collapse = ""),
                            mixed = sapply(status_reg33_AmyloGram, function(single_seq)
                              length(unique(single_seq)) > 1
                            ),
                            AmyLoad_et = ets[reg33_AmyloGram],
                            reg33_et = sapply(status_reg33_AmyloGram, decide_reg)
)

table(reg33_al_comp[, c("AmyLoad_et", "reg33_et")])

# write.csv2(data.frame(name = all_names, reg33 = sapply(seq_list, paste0, collapse = "")), 
#            file = "reg33vsAmyLoad.csv",
#            row.names = FALSE, quote = FALSE)
# 
# write.fasta(lapply(seq_list, paste0, collapse = ""), all_names, "reg33.fasta")
compr <- read.csv2("./data/reg33vsAmyLoad.csv", row.names = NULL, stringsAsFactors = FALSE)

