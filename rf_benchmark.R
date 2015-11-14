library(seqinr)
#library(randomForest)
#library(ranger)
#library(bigrf) #windows only
library(biogram)
library(microbenchmark)
library(pbapply)

create_data <- function(n_seq = 500, n = 2) {
  n_seq <- 500
  seqs <- rbind(t(replicate(n_seq, sample(a()[-1], 6, prob = c(rep(0.7, 4), rep(0.1, 16))))),
                t(replicate(n_seq, sample(a()[-1], 6))))
  dat <- data.frame(as.matrix(count_ngrams(seqs, n, a()[-1])), tar = c(rep("pos", n_seq), rep("neg", n_seq)))
}

do_benchmark <- function(dat)
  microbenchmark(randomForest::randomForest(tar ~ ., dat),
                 party::cforest(tar ~ ., dat),
                 ranger::ranger(tar ~ ., dat)), 
                 #bigrf::bigrfc(dat[, -ncol(dat)], dat[, ncol(dat)]),
                 unit = "s")

nn <- expand.grid(n_seq = c(100, 500, 1000, 2000, 5000), n = 1L:4)

rf_benchmark <- pblapply(1L:nrow(nn), function(i) {
  do_benchmark(create_data(nn[i, 1], nn[i, 2]))
})

save(rf_benchmark, file = "rf_benchmark.RData")
