#read reg33 file

library(dplyr)
library(seqinr)


r33_raw <- readLines("./data/reg33.txt") %>% matrix(ncol = 3, byrow = TRUE) %>% 
  data.frame(stringsAsFactors = FALSE)
amyld_ids <- lapply(strsplit(r33_raw[, 3], " "), function(i) matrix(as.numeric(i), nrow = 2))
#check the lengths
whole_seqs <- strsplit(r33_raw[[2]], "")

seq_list <- unlist(lapply(1L:length(amyld_ids), function(i) {
  lapply(1L:ncol(amyld_ids[[i]]), function(j) {
    id_row <- amyld_ids[[i]][, j]
    whole_seqs[[i]][id_row[1]:id_row[2] + 1]
  })
}), recursive = FALSE)

all_names <- unlist(lapply(1L:length(r33_raw[[1]]), function(i) {
  names_list <- rep(r33_raw[[1]][i], ncol(amyld_ids[[i]]))
  if(ncol(amyld_ids[[i]]) > 1) {
    paste0(names_list, "_", 1L:ncol(amyld_ids[[i]]))
  } else {
    names_list
  }
}))

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]


# write.csv2(data.frame(name = all_names, rep33 = sapply(seq_list, paste0, collapse = "")), 
#            file = "rep33vsAmyLoad.csv",
#            row.names = FALSE, quote = FALSE)
# 
# write.fasta(lapply(seq_list, paste0, collapse = ""), all_names, "reg33.fasta")
compr <- read.csv2("./data/rep33vsAmyLoad.csv", row.names = NULL, stringsAsFactors = FALSE)

