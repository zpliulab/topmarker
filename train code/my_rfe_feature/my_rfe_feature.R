setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_rfe_feature')
library(caret)
library(kernlab)
library(randomForest)

#使用rfe函数进行特征选择

DataPath <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_rfe_feature\\Cluster_expr_label"
FileNames <-dir(DataPath)

for(i in 1:(length(FileNames))){
  data <- read.table(file = paste(FileNames[i],sep="\t"),header = T,sep='\t')
 
  n<-dim(data)[2]
  set.seed(5)
  rfe_comm<-rfe(x=data[,-n],y=as.factor(data[,n]),sizes = c(1:40),
                rfeControl = rfeControl(functions = rfFuncs, method = "cv"))
  
  # p<-plot(rfe_comm, type = c('g','o'),main=paste('Accuracy of ',FileNames[i],sep=""))
  # 
  # pdf(file=paste('Accuracy of ',FileNames[i],'.pdf',sep=''))
  # plot(p)
  # dev.off()
  
  rfe_data<-data[,rfe_comm$optVariables]
  write.table(rfe_data,file=paste('SelectbyRFEofFeatureData_',FileNames[i],sep=""),sep = '\t',quote=F)
  rfe_gene<-data.frame(rfe_comm$optVariables)
  colnames(rfe_gene)<-"gene_select"
  write.table(rfe_gene,file=paste('Genes_SelectbyRFEofFeatureData_',FileNames[i],sep=""),sep = '\t',quote=F,row.names = F)
}


#整理一下筛选出的基因，合并到一个文件中。20个社群共筛选出220基因。
GenePath <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_rfe_feature\\GenesSelect_byRFE"
GeneFileName <- dir(GenePath)
Num<-c()
rfe_gene<-data.frame()
for(k in 1:(length(GeneFileName))){
  data <- read.table(file = paste(GeneFileName[k],sep="\t"),header = T)
  t<-dim(data)[1]
  Num<-c(Num,t)
  clu<-rep(GeneFileName[k],t)
  oneclu<-cbind(data,clu)
  colnames(oneclu)<-c('gene','cluster')
  rfe_gene<-rbind(rfe_gene,oneclu)
}  
write.table(rfe_gene,"rfe_feature_gene.txt",quote=F,sep="\t",row.names = F)

gene_num<-sum(Num)
gene_num
Numf<-data.frame(Num)
GeneFileNamef<-data.frame(GeneFileName)
Num_Clu_RFE<-cbind(GeneFileNamef,Numf)
colnames(Num_Clu_RFE)<-c("Cluster","RFEGene_Num")
write.table(Num_Clu_RFE,"Num_Clu_RFE_feature.txt",quote=F,sep="\t")


#####
#看rfe筛选出的基因与kegg_hcc_pathway中的基因有多少重合
library(tidyr)
rfe_gene<-read.table("rfe_feature_gene.txt",header = TRUE,sep="\t")
rg<-as.data.frame(rfe_gene[,1])
colnames(rg)<-"gene"

hcc_pathway<-read.table("kegg_hcc_pathway.txt",header = TRUE,sep="\t")
colnames(hcc_pathway)<-c("ID","gene")
hp<-separate(hcc_pathway,col="gene",into=c("gene","v1"),sep = ";",remove=T,fill="right")
hp_gene<-data.frame(hp[,2])
colnames(hp_gene)<-"hp_gene"

share_gene<-merge(rg,hp_gene,by.x="gene",by.y="hp_gene",all=FALSE)
dim(share_gene)
share_gene
#167个kegg—hcc-gene，220个差异基因，只有2个交叉···share_gene:AKT2,PRKCG


#33个hcc相关基因与rfe筛选出的基因有没有交叉。
homo_sapiens<-read.table(file="Homo_sapiens.gene_info",sep="\t",header=T,fill=TRUE,quote="")
homo<-subset.data.frame(homo_sapiens,select=c("GeneID","Symbol"))
rfe_gene_ID<-merge(rfe_gene,homo,by.x="gene",by.y="Symbol",all=FALSE)

disease33<-read.table("disease_33.txt",header=F,sep="\t")
dis_33gene<-data.frame(disease33$V2)
colnames(dis_33gene)<-"dis_geneid"
over_gene<-merge(dis_33gene,rfe_gene_ID,by.x="dis_geneid",by.y="GeneID",all=FALSE)
#write.table(over_gene,"RFEgene_33HCCgene_overlap.txt",quote=F,sep="\t",row.names = F)#no overlap gene





#####
#看筛选出的基因pvalue怎么样，在每个社群中的pvalue的排名怎么样。

genename_pvalue_fdr<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t")
#BigClu<-read.table("BigClu.txt",header = TRUE,sep="\t")
rfe_gene<-read.table("rfe_feature_gene.txt",header = TRUE,sep="\t")

# BigClu_pvalue<-merge(BigClu,genename_pvalue_fdr,by.x="gene",by.y="genename",all.x=TRUE)
# BigClu_pvalue_clu<-BigClu_pvalue[order(BigClu_pvalue$cluster1),]
# BigClu_pvalue_order<-BigClu_pvalue[order(BigClu_pvalue$p_value,decreasing = F),]
# 
# write.table(BigClu_pvalue_clu,"BigClu_pvalue_clu.txt",quote=F,sep="\t",row.names = F)
# write.table(BigClu_pvalue_order,"BigClu_pvalue_order.txt",quote=F,sep="\t",row.names = F)

rfe_gene_pvalue<-merge(rfe_gene,genename_pvalue_fdr,by.x="gene",by.y="genename",all.x=TRUE)
rfe_gene_pvalue_clu<-rfe_gene_pvalue[order(rfe_gene_pvalue$cluster),]
rfe_gene_pvalue_order<-rfe_gene_pvalue[order(rfe_gene_pvalue$p_value,decreasing = F),]

write.table(rfe_gene_pvalue_clu,"rfe_feature_gene_pvalue_clu.txt",quote=F,sep="\t",row.names = F)
write.table(rfe_gene_pvalue_order,"rfe_feature_gene_pvalue_order.txt",quote=F,sep="\t",row.names = F)

#BigClu_pvalue$cluster1==2 & BigClu_pvalue$p_value<0.05 : 67/83
#BigClu_pvalue$p_value<0.05 :702/880
#rfe_gene_pvalue$p_value<0.05 : 
#rfe_gene_pvalue$p_value<0.05 & rfe_gene_pvalue$cluster==2 : 

BigClu2 <- BigClu_pvalue[( BigClu_pvalue$p_value<0.05),]











