source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)


ftraits <- choose_properties()

#time-consuming step commented out. Instead just load the results.
aa_groups <- create_encodings(ftraits)
save(aa_groups, file = "aa_groups.RData")
load("aa_groups.RData")
aa_groups <- vector2list(aa_groups)
#need to add seq

pos_gl <- extract_ngrams(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"), aa_groups)
neg_gl <- extract_ngrams(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"), aa_groups)
