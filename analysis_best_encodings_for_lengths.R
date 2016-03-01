summary(amyloids)

amyloids %>% group_by(len_range) %>%
  summarise(which.max(AUC_mean))

amyloids %>% group_by(len_range) %>%
  summarise(max(AUC_mean))

best_enc <- c((amyloids[amyloids$len_range=="(10,15]",])[135191,"enc"],
(amyloids[amyloids$len_range=="(15,25]",])[27378,"enc"],
(amyloids[amyloids$len_range=="[5,6]",])[151207,"enc"],
(amyloids[amyloids$len_range=="(6,10]",])[113669,"enc"])

distances <- matrix(ncol=4, nrow=4)
for(i in 1:4){
  for(j in 1:4){
    distances[i,j] <- calc_ed(aa_groups[[best_enc[i]]], aa_groups[[best_enc[j]]])
  }
}

#distances are quite big
