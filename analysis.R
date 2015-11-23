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
aa_groups <- vector2list(aa_groups)

seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
               read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
seqs_list <- seqs_list[lengths(seqs_list) > 5]

seqs_m <- t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i)))))

gl <- extract_ngrams(seqs_m, aa_groups)
