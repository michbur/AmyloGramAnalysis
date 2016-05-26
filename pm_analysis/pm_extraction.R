library(seqinr)
library(biogram)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i)))))

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

gl <- lapply(1L:nrow(seqs_m), function(i) {
  do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seqs_m[i, ][!is.na(seqs_m[i, ])], 6, a()[-1])), ""))
})

dat <- do.call(rbind, lapply(gl, function(i) {
  res <- data.frame(i)
  colnames(res) <- paste0("P", 1L:6)
  cbind(n_hex = 1L:nrow(res), len = nrow(res) + 5, res)
}))

final_dat <- cbind(prot_name = unlist(sapply(1L:length(gl), function(i) rep(rownames(seqs_m)[i], nrow(gl[[i]])))),
                   prot_id = unlist(sapply(1L:length(gl), function(i) rep(i, nrow(gl[[i]])))),
                   amyloid =  unlist(sapply(1L:length(gl), function(i) rep(ets[i] == 1, nrow(gl[[i]])))),
                   dat
)

write.csv2(final_dat, file = "./pm_analysis/pm_dat.csv", row.names = FALSE)

write.csv2(final_dat[final_dat[["amyloid"]], -3], file = "./pm_analysis/pm_dat_pos.csv", row.names = FALSE)
write.csv2(final_dat[!final_dat[["amyloid"]], -3], file = "./pm_analysis/pm_dat_neg.csv", row.names = FALSE)
write.csv2(aaprop[c("ARGP820101", "BHAR880101", "CHAM820101", "KIMC930101"), ], file = "./pm_analysis/aaprop.csv")
