setwd('D:\\A_My_Data\\HCC\\1_data\\kegg_33gene')
library(igraph)
library(tidyr)

gene33<-read.table("Genes_SelectbyRFEof_expr_label_cluster 3 .txt",header = TRUE,sep="\t")
hcc_ppi<-read.table("hcc_genepairs.txt",header = TRUE,sep="\t")


hcc_pathway<-read.table("kegg_hcc_pathway.txt",header = TRUE,sep="\t")
colnames(hcc_pathway)<-c("ID","gene")
hp<-separate(hcc_pathway,col="gene",into=c("gene","v1"),sep = ";",remove=T,fill="right")
hp_gene<-data.frame(hp[,2])
colnames(hp_gene)<-"gene_select"
write.table(hp_gene,"kegg_gene.txt",quote=F,sep="\t",row.names = F)

kegg_33gene<-rbind(gene33,hp_gene)

share_ppi_1<-merge(kegg_33gene,hcc_ppi,by.x="gene_select",by.y="Gene1",all=FALSE)
share_ppi_2<-merge(kegg_33gene,hcc_ppi,by.x="gene_select",by.y="Gene2",all=FALSE)

colnames(share_ppi_2)<-colnames(share_ppi_1)
share_ppi<-rbind(share_ppi_1,share_ppi_2)

write.table(share_ppi,"share_ppi.txt",quote=F,sep="\t",row.names = F)


kegg_2<-cbind(hp_gene,1)
gene33_2<-cbind(gene33,2)
colnames(kegg_2)<-colnames(gene33_2)

kegg_33gene_2<-rbind(kegg_2,gene33_2)
write.table(kegg_33gene_2,"kegg_33gene_2.txt",quote=F,sep="\t",row.names = F)



#33个hcc相关基因与社区基因有没有交叉。

disease33<-read.table("disease_33.txt",header=F,sep="\t")
homo<-read.table(file="homo.txt",sep="\t",header=T)

dis_33_id<-data.frame(disease33$V2)
colnames(dis_33_id)<-"gene_select_id"
dis_gene<-merge(dis_33_id,homo,by.x="gene_select_id",by.y="GeneID",all=FALSE)
dis_symbol<-data.frame(dis_gene[,2])
colnames(dis_symbol)<-"gene_select"
dis_33gene<-rbind(gene33,dis_symbol)

write.table(dis_symbol,"dis_symbol.txt",quote=F,sep="\t",row.names = F,col.names = F)


dis_ppi_1<-merge(dis_33gene,hcc_ppi,by.x="gene_select",by.y="Gene1",all=FALSE)
dis_ppi_2<-merge(dis_33gene,hcc_ppi,by.x="gene_select",by.y="Gene2",all=FALSE)

colnames(dis_ppi_2)<-colnames(dis_ppi_1)
dis_ppi<-rbind(dis_ppi_1,dis_ppi_2)

write.table(dis_ppi,"dis_ppi.txt",quote=F,sep="\t",row.names = F)


dis_2<-cbind(dis_symbol,1)
gene33_2<-cbind(gene33,2)
colnames(dis_2)<-colnames(gene33_2)

dis_33gene_2<-rbind(dis_2,gene33_2)
write.table(dis_33gene_2,"dis_33gene_2.txt",quote=F,sep="\t",row.names = F)

#####
adj1_node<-read.csv("adj1_node.csv",header=T,sep=",")
adj1_node_2<-data.frame(adj1_node[,3])
write.table(adj1_node_2,"adj1_node.txt",quote=F,sep="\t",row.names = F,col.names = F)

adj1_node_3<-cbind(adj1_node_2,3)
colnames(adj1_node_3)<-colnames(gene33_2)
adj1_node_clu<-rbind(dis_33gene_2,adj1_node_3)
write.table(adj1_node_clu,"adj1_node_clu.txt",quote=F,sep="\t",row.names = F)


dis_ppi_node<-read.csv("dis_ppi_node.csv",header=T,sep=",")
dis_ppi_node_2<-data.frame(dis_ppi_node[,3])
colnames(dis_ppi_node_2)<-"name"

adj2_1<-merge(dis_ppi_node_2,hcc_ppi,by.x="name",by.y="Gene1",all=FALSE)
adj2_2<-merge(dis_ppi_node_2,hcc_ppi,by.x="name",by.y="Gene2",all=FALSE)

colnames(adj2_2)<-colnames(adj2_1)
dis_ppi<-rbind(adj2_1,adj2_2)

write.table(dis_ppi,"adj2_ppi.txt",quote=F,sep="\t",row.names = F)

hcc_adj1_bio_node<-read.csv("hcc_adj1_bio_node.csv",header=T,sep=",")
all_node<-data.frame(hcc_adj1_bio_node$name)
colnames(all_node)<-"node"
write.table(all_node,"all_adj1_node.txt",quote=F,sep="\t",row.names = F)














