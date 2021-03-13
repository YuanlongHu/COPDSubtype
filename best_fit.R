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


Best_fit <- function(res){
  best_k <- cola::suggest_best_k(res)
  if(nrow(best_k[best_k$`1-PAC`>0.9,])==0){
    best_k2 <- apply(best_k[,-c(1,5)],2, function(x){
    ifelse(x == max(x),1,0)
   }) %>%
    as.data.frame()
   best_k2$vote <- apply(best_k2,1,sum)
   best_k2 <- best_k2[best_k2$vote %in% max(best_k2$vote),]
   best_k <- best_k[rownames(best_k2),]
   
   message("Fulfill third rule:\n")
   message(paste(">>> Best-fit:", rownames(best_k2)))
   message(">>> Parameters: ")
   message(paste0("  beat k: ", best_k[1,1], "\n",
                  "  1-PAC: ", best_k[1,2], "\n",
                  "  mean silhouette: ", best_k[1,3], "\n",
                  "  concordance: ", best_k[1,4]))
   message(">>> Vote:")
   message(paste0("  1-PAC: ", best_k2[1,1], "\n",
                  "  mean silhouette: ", best_k2[1,2], "\n",
                  "  concordance: ", best_k2[1,3]))
   }else{
    best_k <- best_k[best_k$`1-PAC`>0.9,]
    message("Fulfill second rule")
    return(best_k)
  }

}


Best_fit(rl)



