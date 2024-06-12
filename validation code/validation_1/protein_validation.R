setwd('D:\\A_My_Data\\HCC\\1_data\\protein_validation')
library(xlsx)

protein<-read.xlsx("nature26140-s6.xlsx",sheetName = "Sheet1", encoding = 'UTF-8')
homo<-read.table("homo.txt",header=T,sep="\t")


protein_clu3<-merge(protein,clu3,by.x="Gene.Name",by.y="gene_select",all=FALSE)


upper_protein<-data.frame(toupper(protein$Gene.Name))
colnames(upper_protein)<-"genename"
protein_clu3<-merge(upper_protein,clu3,by.x="genename",by.y="gene_select",all=FALSE)

protein_homo<-merge(upper_protein,homo,by.x="genename",by.y="Symbol",all=FALSE)


all_adj1_node<-read.table("all_adj1_node.txt",header=T,sep="\t")
protein_all_adj1_node<-merge(upper_protein,all_adj1_node,by.x="genename",by.y="node",all=FALSE)





mouse<-read.table("mouse.txt",header=T,sep="\t",quote = "")#,stringsAsFactors=FALSE
protein_id<-merge(protein,mouse,by.x="Gene.Name",by.y="Symbol",all=FALSE)
protein_homo<-merge(protein_id,homo,by.x="GeneID",by.y="GeneID",all=FALSE)



