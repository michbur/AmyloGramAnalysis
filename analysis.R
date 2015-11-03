source("./functions/choose_properties.R")
source("./functions/create_encodings.R")

ftraits <- choode_properties()
aa_groups <- create_encodings(ftraits)
save(aa_groups, file = "aa_groups.RData")