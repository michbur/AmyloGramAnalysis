source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")

ftraits <- choose_properties()

#time-consuming step commented out. Instead just load the results.
#aa_groups <- create_encodings(ftraits)
#save(aa_groups, file = "aa_groups.RData")
load("aa_groups.RData")
aa_groups <- vector2list(aa_groups)
#need to add seq
extract_ngrams(seq, aa_groups)
