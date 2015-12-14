# report preamble

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(cvTools)
require(ranger)

# functions

format_property_table <- function(ftraits) {
  add_names <- c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
                 "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")
  traits_names <- c(sapply(aaindex, function(i) i[["D"]]), 
                    add_names)
  unname(traits_names[ftraits])
  
  prop_id <- c(2, 5, 7, 12, 16, 17, #contactivity,
               4, 6, 10, #beta freq
               1, 9, 14, 15, #hyrophobicity
               3, 8, #polarity
               13, #size
               11) #stability
  
  # sequence as given by prop_id
  prop_category <- c(rep("Contactivity", 6), 
                     rep("$\\beta$-frequency", 3),
                     rep("Hydrophobicity", 4), 
                     rep("Polarity", 2), 
                     "Size",
                     "Stability")
  
  data.frame(Category = prop_category, Property = unname(traits_names[ftraits])[prop_id])
}