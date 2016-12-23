library(Biostrings)

remove_dash <- function(x)
  x[x != "-"]

seq2bio <- function(x)
  AAString(paste0(x, collapse = ""))

all_prots <- c(read.fasta("/data/amyloid_neg_benchmark.fasta", seqtype = "AA"), 
               read.fasta("/data/amyloid_pos_benchmark.fasta", seqtype = "AA")) %>% 
  lapply(toupper) %>% 
  lapply(remove_dash) %>% 
  lapply(seq2bio)

pep424 <- read.fasta("/benchmark/pep424_better_names.fasta", seqtype = "AA") %>% 
  lapply(seq2bio)

pep424_pid <- lapply(pep424, function(i) 
  pbsapply(all_prots, function(j) 
    pid(pairwiseAlignment(i, j, type="global"), type = "PID1")
  )
)

save(pep424_pid, file = "/results/pep424_pid.RData")
