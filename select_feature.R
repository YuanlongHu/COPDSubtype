devtools::install_github("YuanlongHu/yltool")
library(yltool)
library(Boruta)
library(magrittr)
library(UpSetR)
library(VennDiagram)
#####################################################
# ECLIPSE_data: Gene expression data of ECLIPSE
# ECLIPSE_pdata: Phenotype data of ECLIPSE
res_list_ECLIPSE_FEV1 <- yltool::lm_single(ECLIPSE_data, ECLIPSE_pdata$FEV1)
res_list_ECLIPSE_FEV1FVC <- yltool::lm_single(ECLIPSE_data, ECLIPSE_pdata$FEV1FVC)
res_list_ECLIPSE_FEV1$p.adjust <- p.adjust(res_list_ECLIPSE_FEV1$Pvalue, method = "BH")
res_list_ECLIPSE_FEV1FVC$p.adjust <- p.adjust(res_list_ECLIPSE_FEV1FVC$Pvalue, method = "BH")

res_list_ECLIPSE_FEV1 <- res_list_ECLIPSE_FEV1[res_list_ECLIPSE_FEV1$Pvalue<0.05,]
res_list_ECLIPSE_FEV1FVC <- res_list_ECLIPSE_FEV1FVC[res_list_ECLIPSE_FEV1FVC$Pvalue<0.05,]

# Boruta with "p-value=0.01"
set.seed(12345)
res_Boruta_ECLIPSE_FEV1 <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1$y,], 
                                                 pdata=ECLIPSE_pdata$FEV1,
                                                 pValue = 0.01)
res_Boruta_ECLIPSE_FEV1FVC <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1FVC$y,], 
                                                 pdata=ECLIPSE_pdata$FEV1FVC,
                                                 pValue = 0.01)
gene_s2 <- union(getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC),
                 getSelectedAttributes(res_Boruta_ECLIPSE_FEV1))

########################################################################
# Boruta with "p-value=0.05"
set.seed(12345)
res_Boruta_ECLIPSE_FEV1_t1 <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1$y,], 
                                                 pdata=ECLIPSE_pdata$FEV1,
                                                 pValue = 0.05)
res_Boruta_ECLIPSE_FEV1FVC_t1 <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1FVC$y,], 
                                                    pdata=ECLIPSE_pdata$FEV1FVC,
                                                    pValue = 0.05)
Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC_t1)
Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC)

gene_s_t1 <- union(Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC_t1),
      Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1_t1))

# Boruta with "p-value=0.1"
set.seed(12345)
res_Boruta_ECLIPSE_FEV1_t2 <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1$y,], 
                                                    pdata=ECLIPSE_pdata$FEV1,
                                                    pValue = 0.1)
res_Boruta_ECLIPSE_FEV1FVC_t2 <- yltool::select_Boruta(expr=ECLIPSE_data[res_list_ECLIPSE_FEV1FVC$y,], 
                                                       pdata=ECLIPSE_pdata$FEV1FVC,
                                                       pValue = 0.1)
union(Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC_t2),
      Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1_t2)) %>%
  intersect(gene_s2)

gene_s_t2 <- union(Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1FVC_t2),
      Boruta::getSelectedAttributes(res_Boruta_ECLIPSE_FEV1_t2))

##############################################################

res_list_ECLIPSE_FEV1 <- res_list_ECLIPSE_FEV1[order(abs(res_list_ECLIPSE_FEV1$Beta), decreasing = T),]

res_list_ECLIPSE_FEV1FVC <- res_list_ECLIPSE_FEV1FVC[order(abs(res_list_ECLIPSE_FEV1FVC$Beta), decreasing = T),]
res_list_ECLIPSE_FEV1FVC$y[1:round(nrow(res_list_ECLIPSE_FEV1FVC)*0.1)]

gene_s_t3 <- union(
  res_list_ECLIPSE_FEV1$y[1:round(nrow(res_list_ECLIPSE_FEV1)*0.1)],
  res_list_ECLIPSE_FEV1FVC$y[1:round(nrow(res_list_ECLIPSE_FEV1FVC)*0.1)])
gene_s_t3


gene_s_t4 <- union(
  res_list_ECLIPSE_FEV1$y[1:round(nrow(res_list_ECLIPSE_FEV1)*0.05)],
  res_list_ECLIPSE_FEV1FVC$y[1:round(nrow(res_list_ECLIPSE_FEV1FVC)*0.05)])
gene_s_t4


list(Boruta_0.01 = gene_s2,
     Boruta_0.05 = gene_s_t1,
     Boruta_0.1 = gene_s_t2,
     LM_0.1 = gene_s_t3,
     LM_0.05=gene_s_t4
) %>%
  venn.diagram(
    filename = NULL,col = "black",lty = "dotted",lwd = 2, 
    fill = c("cornflowerblue", "yellow", "darkorchid1","red","white"),
    alpha = 0.50,cex = 1,fontfamily = "serif",fontface = "bold",
    cat.col = c("black", "black", "black","black","black"),
    cat.cex = 1.8,cat.fontface = "bold", cat.fontfamily = "serif",margin = 0.1
  ) %>%
  grid.draw()


upset(fromList(list(Boruta_0.01 = gene_s2,
                    Boruta_0.05 = gene_s_t1,
                    Boruta_0.1 = gene_s_t2,
                    LM_0.1 = gene_s_t3,
                    LM_0.05=gene_s_t4
)), order.by = "freq", empty.intersections = "on")
