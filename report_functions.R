source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")
source("./functions/plot_tools.R")

if(Sys.info()["nodename"] == "tobit" )
  setwd("~/Dropbox/doktorat/moje_prace/amyloid_ngram2/")


library(dplyr)
library(biogram) #calc_ed
require(seqinr) #for choose_properties


# read data ----------------------------------------

amyloids <- read.csv("results/all_summaries.csv")
load("aa_groups.RData")
aa_groups <- string2list(aa_groups)


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

trait_tab <- read.csv(file = "./results/trait_tab.csv")

# best ecodings ----------------------------------------------

best_positions <- amyloids %>% 
  select(len_range, enc_adj, AUC_mean) %>%
  group_by(len_range) %>%
  mutate(position = rank(AUC_mean)) %>%
  ungroup %>%
  group_by(enc_adj) %>%
  summarise(cum_rank = sum(position)) %>%
  arrange(desc(cum_rank)) %>%
  slice(1)

best_enc <- best_positions[["enc_adj"]]

load("./results/cv_results_full.RData")

full_alphabet <- lapply(1L:length(cv_results_full[[1]]), function(single_replicate_id) {
  lapply(1L:5, function(single_fold) 
    data.frame(cv_results_full[[1]][[single_replicate_id]][[single_fold]][[1]], replicate = floor((single_replicate_id - 1)/15) + 1)
  )
}) %>% unlist(recursive = FALSE) %>%
  do.call(rbind, .) %>% 
  mutate(MCC = calc_mcc(TP, TN, FP, FN)) %>%
  select(replicate, len_range, AUC, MCC, Sens, Spec) %>%
  group_by(replicate, len_range) %>%
  summarize_each(funs(mean = liberal_mean, sd = liberal_sd), AUC, MCC, Sens, Spec) %>%
  ungroup %>% 
  mutate(pos = c(6, 10, 15)[replicate])


# Fig 1 all encodings sens/spec  ----------------------------------------

amyloids_plot <- select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
  mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
  select(AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, et, enc_adj) %>%
  rbind(select(full_alphabet, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>% 
          mutate(et = "full alphabet", enc_adj = 0)) %>%
  mutate(len_range = factor(len_range, levels = c("[5,6]", "(6,10]", "(10,15]", "(15,25]")),
         pos = factor(pos, labels = paste0("Training peptide\nlength: ", 
                                           c("6", "10 or less", "15 or less"))),
         et = factor(et, labels = c("Reduced alphabet", "Best performing reduced alphabet",
                                    "Reduced alphabet from literature", "Full alphabet"))) %>%
  mutate(len_range = factor(len_range, 
                            labels = paste0("Testing peptide\nlength: ", c("6 ", "7-10", "11-15", "16-25"))))

sesp_plot <- ggplot(amyloids_plot, aes(x = Sens_mean, y = Spec_mean, color = et, shape = et)) +
  geom_point() +
  scale_color_manual("", values = c("grey", "red", "blue", "green")) +
  scale_shape_manual("", values = c(1, 16, 15, 15)) +
  scale_x_continuous("Mean sensitivity") +
  scale_y_continuous("Mean specificity") +
  facet_grid(pos ~ len_range) +
  my_theme +
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
             aes(x = Sens_mean, y = Spec_mean, color = et))


#png("./pub_figures/sesp_plot.png", height = 648, width = 648)
cairo_pdf("./pub_figures/sesp_plot.pdf", height = 9, width = 9)
print(sesp_plot)
dev.off()


AUC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = AUC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1) +
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
             aes(x = len_range, y = AUC_mean, color = et, shape = et)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean AUC") +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "red", "blue", "green"), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip()

cairo_pdf("./pub_figures/AUC_boxplot.pdf", height = 9, width = 9)
#png("./pub_figures/AUC_boxplot.png", height = 648, width = 648)
print(AUC_boxplot)
dev.off()

MCC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = MCC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1) +
  geom_point(data = filter(amyloids_plot, et != "Reduced alphabet"), 
             aes(x = len_range, y = MCC_mean, color = et, shape = et)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean MCC") +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "red", "blue", "green"), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip()

cairo_pdf("./pub_figures/MCC_boxplot.pdf", height = 9, width = 9)
#png("./pub_figures/MCC_boxplot.png", height = 648, width = 648)
MCC_boxplot
dev.off()
