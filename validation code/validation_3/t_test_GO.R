setwd('D:\\A_My_Data\\HCC\\1_data\\validation_3')
library(tidyverse)
library(httr)
library(jsonlite)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(pathview)

#data <- read.table(file = "Genes_SelectbyRFEof_expr_label_cluster 3 .txt",sep="\t",header = T)

homo<-read.table(file="homo.txt",sep="\t",header=T)
#gpf2<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t")
gpf<-read.table("p_limma_order.txt",header = TRUE,sep="\t")

data<-data.frame(gpf$genename[1:33])
colnames(data)<-"genename"

gene_id_sym <- merge(data,homo,by.x="genename",by.y="Symbol",all=FALSE)#_select
gene_id<-gene_id_sym[,2]#取出每个社群中包含的基因名，匹配成基因ID，变成逗号相连的list。


ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = gene_id,
                   pvalueCutoff = 0.05,pAdjustMethod = "fdr",
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- as.data.frame(ego_BP)


ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = gene_id,
                   pvalueCutoff = 0.05,pAdjustMethod = "fdr",
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)


ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = gene_id,
                   pvalueCutoff = 0.05,pAdjustMethod = "fdr",
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)

go_result_df<-rbind(ego_result_BP,ego_result_CC,ego_result_MF)

write.table(ego_result_BP,"t_test_ego_BP.txt",quote = F,sep='\t',row.names = F)
write.table(ego_result_CC,"t_test_ego_CC.txt",quote = F,sep='\t',row.names = F)
write.table(ego_result_MF,"t_test_ego_MF.txt",quote = F,sep='\t',row.names = F)
write.table(go_result_df,"t_test_ego_df.txt",quote = F,sep='\t',row.names = F)


# go_enrich_df <- data.frame(       ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
#                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
#                                   GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
#                                   Pvalue=c(ego_result_BP$pvalue, ego_result_CC$pvalue, ego_result_MF$pvalue),
#                                   geneID=c(ego_result_BP$geneID, ego_result_CC$geneID, ego_result_MF$geneID),
#                                   type=factor(c(rep("biological process", dim(ego_result_BP)[1]), 
#                                                 rep("cellular component", dim(ego_result_CC)[1]),
#                                                 rep("molecular function", dim(ego_result_MF)[1])
#                                   )))#levels=c("biological process", "cellular component", "molecular function")
# 
# write.table(go_enrich_df,"rfe_clu3_GO.txt",quote = F,sep='\t',row.names = F)

  
ee2<-emapplot(ego_BP,title="rfe_clu3_GO")
plot(ee2)

data(geneList, package="DOSE")

cnetplot(ego_BP, foldChange=geneList,showCategory = 80)#
heatplot(ego_BP, foldChange=geneList,showCategory = 80)
heatplot(ego_CC, foldChange=geneList,showCategory = 80)
heatplot(ego_MF, foldChange=geneList,showCategory = 80)


heatplot(ego_BP, foldChange=geneList)
heatplot(ego_CC)


#kegg
# data <- read.table(file = "Genes_SelectbyRFEof_expr_label_cluster 3 .txt",sep="\t",header = T)
# gene_id_sym <- merge(data,homo,by.x="gene_select",by.y="Symbol",all=FALSE)
# gene_id<-gene_id_sym[,2]#取出每个社群中包含的基因名，匹配成基因ID，变成逗号相连的list。

kegg <- enrichKEGG(gene = gene_id,
                   organism = 'hsa', #KEGG可以用organism = 'hsa'
                   pvalueCutoff = 1)
#setReadable(kegg, OrgDb = org.Hs.eg.db,keyType="SYMBOL")

keggf<-data.frame(kegg)
write.table(keggf,file="rfe_clu3_KEGG.txt",quote = F,sep='\t',row.names = F)


kk<-dotplot(kegg,title="rfe_clu3_KEGG")
plot(kk)

ee<-emapplot(kegg,title="rfe_clu3_KEGG")
plot(ee)
heatplot(kegg)

