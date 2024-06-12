setwd("C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_clusterProfiler\\my_allCluster_GO_2")
library(tidyverse)
library(httr)
library(jsonlite)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(pathview)
#对rfe之前的20个社群中的基因进行GO分析，此时有880个基因。

homo_sapiens<-read.table(file="Homo_sapiens.gene_info",sep="\t",header=T,fill=TRUE,quote="")
homo<-subset.data.frame(homo_sapiens,select=c("GeneID","Symbol"))

path <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_clusterProfiler\\my_allCluster_GO_2\\ClusterData"
FileNames <- dir(path)
######  
#对每个社群进行GO分析并画图。  
for(i in 1:(length(FileNames))){#
  data <- read.table(file = paste(FileNames[i],sep = "\t"),header = T)
  gene_id_sym <- merge(data,homo,by.x="gene",by.y="Symbol",all=FALSE)
  gene_id<-gene_id_sym[,2]#取出每个社群中包含的基因名，匹配成基因ID，变成逗号相连的list。
  
  # GO enrichment with clusterProfiler
  # ego_all <- enrichGO(OrgDb="org.Hs.eg.db",
  #                    gene = gene_id,
  #                    pvalueCutoff = 0.05,
  #                    ont = "BP",
  #                    readable=TRUE)
  # ego_result_all <- as.data.frame(ego_all)#1:display_number[1],
  
  ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                     gene = gene_id,
                     pvalueCutoff = 0.05,
                     ont = "BP",
                     readable=TRUE)
  ego_result_BP <- as.data.frame(ego_BP)[]#1:display_number[3],
  if(dim(ego_result_BP)[1]>10){  ego_result_BP <- ego_result_BP[1:10,]  }
  if(dim(ego_result_BP)[1]==0){  ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                                                    gene = gene_id,
                                                    pvalueCutoff = 0.2,
                                                    ont = "MF",
                                                    readable=TRUE)
  ego_result_BP <- as.data.frame(ego_BP)[1,]  }
  # ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]
  
  ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                     gene = gene_id,
                     pvalueCutoff = 0.05,
                     ont = "CC",
                     readable=TRUE)
  ego_result_CC <- as.data.frame(ego_CC)[]#1:display_number[2],
  if(dim(ego_result_CC)[1]>10){  ego_result_CC <- ego_result_CC[1:10,]  }
  if(dim(ego_result_CC)[1]==0){  ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                                                    gene = gene_id,
                                                    pvalueCutoff = 0.2,
                                                    ont = "MF",
                                                    readable=TRUE)
  ego_result_CC <- as.data.frame(ego_CC)[1,]  }
  # ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
  
  
  ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                     gene = gene_id,
                     pvalueCutoff = 0.05,
                     ont = "MF",
                     readable=TRUE)
  ego_result_MF <- as.data.frame(ego_MF)[]#1:display_number[1],
  if(dim(ego_result_MF)[1]>10){  ego_result_MF <- ego_result_MF[1:10,]  }
  if(dim(ego_result_MF)[1]==0){  ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                                                    gene = gene_id,
                                                    pvalueCutoff = 0.2,
                                                    ont = "MF",
                                                    readable=TRUE)
  ego_result_MF <- as.data.frame(ego_MF)[1,]  }
  # ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
  
  
  
  go_enrich_df <- data.frame(       ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                    Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                                    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                    Pvalue=c(ego_result_BP$pvalue, ego_result_CC$pvalue, ego_result_MF$pvalue),
                                    geneID=c(ego_result_BP$geneID, ego_result_CC$geneID, ego_result_MF$geneID),
                                    type=factor(c(rep("biological process", dim(ego_result_BP)[1]), 
                                                  rep("cellular component", dim(ego_result_CC)[1]),
                                                  rep("molecular function", dim(ego_result_MF)[1])
                                    )))#levels=c("biological process", "cellular component", "molecular function")
  
  write.table(go_enrich_df,file=paste("GO",FileNames[i],sep=""),quote = F,sep='\t',row.names = F)
  
  ## numbers as data on x axis
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  ## shorten the names of GO terms
  shorten_names <- function(x, n_word=4, n_char=40){
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
    {
      if (nchar(x) > 40) x <- substr(x, 1, 40)
      x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
      return(x)
    }
    else
    {
      return(x)
    }
  }
  
  labels=(sapply(
    levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
    shorten_names))
  names(labels) = rev(1:nrow(go_enrich_df))
  
  ## colors for bar // green, blue, orange
  CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
  
  pv<-data.frame(go_enrich_df$Pvalue,go_enrich_df$number,go_enrich_df$type)
  colnames(pv)<-c('Pvalue','number','type')
  pv$Pvalue<-signif(pv$Pvalue,digits = 3)
  p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber,fill=type)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() + 
    scale_fill_manual(values = CPCOLS) + theme_bw() + 
    scale_x_discrete(labels=labels) +
    geom_text(data=pv,mapping=aes(x=number, y=1),label=pv$Pvalue,size=2.8)+#标注p值
    xlab("GO term") + 
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    labs(title = paste("The Most Enriched GO Terms of",FileNames[i],sep=" "))
  
  pdf(file=paste("GO_enrichment_",FileNames[i],".pdf",sep=""))
  plot(p)
  dev.off()
  
}


#20个社群每个rfe筛选之后剩多少个基因。将这些基因写出。
######
# fileName <- dir(pattern='[0-9] cluster.txt')
# num<-c()
# rfe_gene1<-data.frame()
# for(k in 1:length(fileName)){
#   data <- read.table(file = paste(fileName[k],sep="\t"),header = T)
#   t<-dim(data)[1]
#   num<-c(num,t)
#   rfe_gene1<-rbind(rfe_gene1,data)
# }  
# fileName2 <- dir(pattern='[0-9]{2} cluster.txt')
# num2<-c()
# rfe_gene2<-data.frame()
# for(k in 1:length(fileName2)){
#   data <- read.table(file = paste(fileName2[k],sep="\t"),header = T)
#   t<-dim(data)[1]
#   num2<-c(num2,t)
#   rfe_gene2<-rbind(rfe_gene2,data)
# }    
# gene_num<-sum(num,num2)  
# rfe_gene<-rbind(rfe_gene1,rfe_gene2)
# write.table(rfe_gene,"rfe_gene.txt",quote=F,sep="\t",row.names = F)




#####
#对每个社群进行kegg分析并画图。  
for(j in 1:(length(FileNames))){#
  data <- read.table(file = paste(FileNames[j],sep="\t"),header = T)
  gene_id_sym <- merge(data,homo,by.x="gene",by.y="Symbol",all=FALSE)
  gene_id<-gene_id_sym[,2]#取出每个社群中包含的基因名，匹配成基因ID，变成逗号相连的list。
   
  kegg <- enrichKEGG(gene = gene_id,
                     organism = 'hsa', #KEGG可以用organism = 'hsa'
                     pvalueCutoff = 0.05)
  keggf<-data.frame(kegg)
  write.table(keggf,file=paste("kegg_",FileNames[j],sep=''),quote = F,sep='\t',row.names = F)
  
  
  kk<-dotplot(kegg,title=paste("kegg enrichment of ",FileNames[j],sep=''))
  pdf(file=paste("kegg_",FileNames[j],".pdf",sep=''))
  plot(kk)
  dev.off()
  
  ee<-emapplot(kegg,title=paste("kegg enrichment of ",FileNames[j],sep='', showCategory = 30))
  pdf(file=paste("kegg_ema_",FileNames[j],".pdf",sep=''))
  plot(ee)
  dev.off()
  
}




# hsa <- pathview(gene.data = gene_id,
# 
#                      pathway.id = "hsa05168", #上述结果中的hsa04120通路
# 
#                      species = "hsa"
# 
#                      )

# source("https://bioconductor.org/biocLite.R")
# biocLite("DOSE")
# biocLite("topGO")
# biocLite("pathview")



