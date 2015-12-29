# Simulate amyloid sequences

# Assumptions:
# 1. Length of amyloidogenic region: 6 aa
# 2. Only four traits responsible for amyloidogenicity

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(cvTools)
require(ranger)
require(hmeasure)

ftraits <- choose_properties()

# later add sample here
imp_prop <- ftraits[1L:4]

#simulate sequence of length n
sim_seq <- function(n) {
  sample(a()[-1], n, replace = TRUE)
}