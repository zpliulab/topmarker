setwd("C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_community\\my_comm_2")
library(tidyverse)
library(httr)
library(jsonlite)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(topGO)
library(pathview)
#对20个社群进行rfe之后的基因进行GO分析，此时有472个基因。

homo_sapiens<-read.table(file="Homo_sapiens.gene_info",sep="\t",header=T,fill=TRUE,quote="")
homo<-subset.data.frame(homo_sapiens,select=c("GeneID","Symbol"))


######  
#对每个社群进行GO分析并画图。  
#for(i in 1:(length(RFE_FileNames))){#
data <- read.csv("Cluster2node_MCODE.csv",sep = "\t",header = T)
colnames(data)<-c("v1")
genep<-separate(data,col="v1",into=c("a1","a2","a3","a4","name","a6","a7"),sep = ",",remove=T,fill="right")

Gene <- subset.data.frame(genep,select=c("name"))
write.table(Gene,"MCODE_cluster2.txt",quote = F,sep="\t")

gene_id_sym <- merge(Gene,homo,by.x="name",by.y="Symbol",all=FALSE)
gene_id<-gene_id_sym[,2]#取出每个社群中包含的基因名，匹配成基因ID，变成逗号相连的list。


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

"cluster mcode"
write.table(go_enrich_df,file=paste("cluster2 mcode.txt"),quote = F,sep='\t',row.names = F)

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
  labs(title = paste("cluster2 mcode"))


plot(p)
  
  
  
#}