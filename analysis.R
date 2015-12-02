source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)


ftraits <- choose_properties()

#time-consuming step commented out. Instead just load the results.
#aa_groups <- create_encodings(ftraits)
#save(aa_groups, file = "aa_groups.RData")
load("aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
               read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
seqs_list <- raw_seqs_list[lengths(raw_seqs_list) > 5]

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[lengths(raw_seqs_list) > 5]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))

#for single aa_group around 85 sec
extracted_ngrams <- extract_ngrams(seqs_m, aa_groups)
#extracted_ngrams1_2500 <- extract_ngrams(seqs_m, aa_groups[1L:2500])
#save(extracted_ngrams1_2500, file = "extracted_ngrams1_2500.RData")
#cat("Dobra robota", file = "/home/michal/Dropbox/done.txt")
