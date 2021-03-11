library(cola)
library(magrittr)
#######################################################
# fit all consensus clustering models
# ECLIPSE_data: whole ECLIPSE gene expression
# gene_s2: 19 selected genes
rl <- run_all_consensus_partition_methods(
  as.matrix(ECLIPSE_data[gene_s2,]),
  top_value_method = c("SD","CV","MAD","ATC"),
  partition_method = c("hclust","kmeans","pam")
  )

best_k <- cola::suggest_best_k(rl)

best_k[best_k$`1-PAC`>0.9,] # fulfill second rule

best_k <- apply(best_k[,-c(1,5)],2, function(x){  # Third rule
  ifelse(x == max(x),1,0)
}) %>%
  as.data.frame()

best_k$vote <- apply(best_k,1,sum)
best_k[best_k$vote %in% max(best_k$vote),] # fulfill third rule

