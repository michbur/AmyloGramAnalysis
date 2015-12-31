-# Simulate amyloid sequences

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
#u - alphabet
#n - length
sim_seq <- function(u, n) {
  sample(u, n, replace = TRUE)
}

#simulate positive
#ngram - a VECTOR of ngram
sim_pos <- function(u, n, ngram) {
  seq <- sim_seq(u, n)
  ngram_pos <- sample(1L:(length(seq) - length(ngram) + 1), 1)
  seq[ngram_pos:(ngram_pos + length(ngram) -1)] <- ngram
  seq
}


ngrams <- list(c("3", "3", "1"),
               c("3", "3"),
               c("3", "1", "3"),
               c("3", "2", "3"),
               c("3", "3", "3"))

simulated_seqs <- rbind(do.call(rbind, lapply(1L:200, function(dummy) 
  sim_pos(as.character(1L:4), 6, sample(ngrams, 1)[[1]]))),
  do.call(rbind, lapply(1L:200, function(dummy) 
    sim_seq(as.character(1L:4), 6))))
  

bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                         ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                         seq = simulated_seqs,
                                         u = 1L:4))

bitrigrams <- bitrigrams > 0
storage.mode(bitrigrams) <- "integer"

test_bis <- test_features(c(rep(1, 200), rep(0, 200)), bitrigrams, adjust = NULL)
imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
