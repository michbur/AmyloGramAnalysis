source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")


require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)
require(hmeasure)
require(pbapply)

load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 16
seqs_list <- raw_seqs_list[purified_seqs_id]

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

load("./results/alns.RData")

perm_alns <- do.call(cbind, unlist(lapply(c("./results/perm_alns5.txt",
                        "./results/perm_alns5_2.txt",
                        "./results/perm_alns5_3.txt"), function(i)
  lapply(strsplit(readLines(i), ";"), as.numeric)), recursive = FALSE))

p_scores <- unlist(lapply(1L:length(alns), function(i)
  mean(perm_alns[i, ] >= alns[i])))

table(p_scores)
mean(p_scores > 0.8)

# what do we align? 
# 0:non-amyloid/non-amyloid
# 1: amyloid/non-amyloid
# 2: amyloid/amyloid, 

aln_type <- unlist(lapply(combn(1L:length(seqs_list), 2, simplify = FALSE), 
              function(single_combn) sum(ets[single_combn])))

p_score_df <- do.call(rbind, lapply(0L:2, function(single_type)
  cbind(type = single_type,
        as.data.frame(table(p_score = p_scores[aln_type == single_type])))
)) %>% 
  mutate(type = factor(type, labels = c("Non-amyloid/non-amyloid", 
                                        "Amyloid/non-amyloid", 
                                        "Amyloid/amyloid"))) %>% 
  group_by(type) %>% 
  mutate(freq = Freq/sum(Freq)) %>% 
  ungroup

library(ggplot2)

ggplot(p_score_df, aes(x = type, y = freq, fill = p_score, 
                       label = formatC(freq, 4, format = "f"))) +
  geom_bar(stat = "identity") +
  geom_text(position = "stack") +
  scale_y_continuous("Frequency") +
  scale_fill_discrete("P-value") +
  scale_x_discrete("Labels of sequences in the alignment")

