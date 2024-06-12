setwd("D:\\A_My_Data\\HCC\\1_data\\validation_3")


limma_p_fc<-read.table("tumor_control_limma.txt",header=T,sep="\t")
clu3<-read.table("Genes_SelectbyRFEof_expr_label_cluster 3 .txt",header = T,sep = "\t")

l1<-limma_p_fc[abs(limma_p_fc$logFC)>1 & limma_p_fc$P.Value<0.05,]

l11<-limma_p_fc[ limma_p_fc$P.Value<0.01,]
l111<-l11[order(l11$P.Value,decreasing = F),]
l111$genename<-as.character(row.names(l111)) 
l111<-l111[,c(7,1:6)]

write.table(l111,"p_limma_order.txt",quote=FALSE,sep="\t")




l1_gene<-data.frame(row.names(l1))
colnames(l1_gene)<-"gene"

same_gene<-merge(clu3,l1_gene,by.x="gene_select",by.y="gene",all=FALSE)


l2<-cbind(row.names(limma_p_fc),limma_p_fc)
colnames(l2)[1]<-"gene"
gene33_fc_p<-merge(clu3,l2,by.x="gene_select",by.y="gene",all=FALSE)

write.table(gene33_fc_p,"gene33_fc_p_limma.txt",quote=FALSE,sep="\t")

l3<-limma_p_fc[limma_p_fc$P.Value<0.01,]
l3<-l3[order(l3$P.Value,decreasing = F),]
write.table(l3,"gene33_p_limma.txt",quote=FALSE,sep="\t")




