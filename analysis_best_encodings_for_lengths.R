
if(Sys.info()["nodename"] == "tobit" )
  setwd("~/Dropbox/doktorat/moje_prace/amyloid_ngram2/")

amyloids <- read.csv("results/all_summaries.csv")
load("aa_groups.RData")
aa_groups <- string2list(aa_groups)

library(dplyr)
library(biogram) #calc_ed
require(seqinr) #for choose_properties

summary(amyloids)

amyloids %>% group_by(len_range) %>%
  summarise(which.max(AUC_mean))

amyloids %>% group_by(len_range) %>%
  summarise(max(AUC_mean))

best_enc <- c((amyloids[amyloids$len_range=="(10,15]",])[135191,"enc_adj"],
              (amyloids[amyloids$len_range=="(15,25]",])[27378,"enc_adj"],
              (amyloids[amyloids$len_range=="[5,6]",])[151207,"enc_adj"],
              (amyloids[amyloids$len_range=="(6,10]",])[113669,"enc_adj"])

distances <- matrix(ncol=4, nrow=4)
for(i in 1:4){
  for(j in 1:4){
    distances[i,j] <- calc_ed(aa_groups[[best_enc[i]]], aa_groups[[best_enc[j]]])
  }
}

#distances are quite big

load("./results/enc_dupes.RData")
#remember, indices in enc_adj and are equal to indices in enc_dupes + 2 and trait_tab  + 2

all_encs <- unlist(enc_dupes[["aa_duplicates"]][best_enc], use.names = FALSE) %>%
  sapply(function(i) substr(i, 3, nchar(i))) %>%
  strsplit("K") %>%
  lapply(as.numeric) %>%
  do.call(rbind, .) %>%
  data.frame %>%
  select(X1) %>%
  unlist(use.names = FALSE) - 2

trait_tab <-  read.csv(file = "./results/trait_tab.csv")
frequencies <- data.frame(trait_tab[all_encs,] %>%
                            apply(1, na.omit) %>% unlist %>% table)
colnames(frequencies) <- c("X", "Freq")
frequencies$X <- as.numeric(as.character(frequencies$X))

#############this part is copied from function choose_properties ################
prop_MK <- read.csv2("AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1]) %>%
  as.numeric

prop_MK <- cbind(prop_MK, years = years) %>% filter(years >= 1980) %>%
  select(X, name) %>% rbind(., 
                            data.frame(X = 545:550, 
                                       name =  c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")))
#################################################################################
prop_MK %>%
  inner_join(frequencies) %>%
  arrange(desc(Freq)) %>%
  select(name, Freq)

