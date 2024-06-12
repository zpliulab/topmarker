setwd('D:\\A_My_Data\\HCC\\1_data\\validation_3')
library(dplyr)
library(e1071)
library(ROCR)
library(pROC)
library(ggplot2)
library(caret)
library(tidyr)
library(tidyverse)

##t-test order
#genename_pvalue_fdr<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t",check.names=FALSE)
genename_pvalue_fdr<-read.table("p_limma_order.txt",header = TRUE,sep="\t",check.names=FALSE)


##GSE25097
GSE25097_expr <- read.table("GSE25097_expr.txt",sep="\t",header = T)
col_GSE25097_expr<-data.frame(row.names(GSE25097_expr))
colnames(col_GSE25097_expr)<-"gene"
GSE25097_expr2<-cbind(col_GSE25097_expr,GSE25097_expr)


##GSE14520、GSE22058、GSE63898
GSE14520_expr <- read.table("GSE14520_expr.txt",sep="\t",header = T)
col_GSE14520_expr<-data.frame(GSE14520_expr$Symbol)
colnames(col_GSE14520_expr)<-"gene"
gse14520_GSM_label<-read.table("gse14520_GSM_label.txt",header=T,sep="\t")
row.names(gse14520_GSM_label)<-gse14520_GSM_label[,1]

GSE22058_expr <- read.table("GSE22058_expr.txt",sep="\t",header = T)
col_GSE22058_expr<-data.frame(GSE22058_expr$Symbol)
colnames(col_GSE22058_expr)<-"gene"
gse22058_GSM_label<-read.table("gse22058_GSM_label.txt",header=T,sep="\t")
row.names(gse22058_GSM_label)<-gse22058_GSM_label[,1]


GSE63898_expr <- read.table("GSE63898_expr.txt",sep="\t",header = T)
col_GSE63898_expr<-data.frame(GSE63898_expr$Symbol)
colnames(col_GSE63898_expr)<-"gene"
gse63898_GSM_label<-read.table("gse63898_GSM_label.txt",header=T,sep="\t")
row.names(gse63898_GSM_label)<-gse63898_GSM_label[,1]


##GSE64041、GSE45436、RNA_seq
GSE64041_expr <- read.table("GSE64041_expr.txt",sep="\t",header = T)
col_GSE64041_expr<-data.frame(GSE64041_expr$symbol)
colnames(col_GSE64041_expr)<-"gene"

GSE45436_expr <- read.table("GSE45436_expr.txt",sep="\t",header = T)
col_GSE45436_expr<-data.frame(GSE45436_expr$Symbol)
colnames(col_GSE45436_expr)<-"gene"
GSE45436_mylabel<-read.table("GSE45436_mylabel.txt",header=T,sep="\t")


LIHC_gene2 <- read.table("LIHC_2_gene.normalized_RNAseq__tissueTypeAll__20180815160255.txt",sep="\t",
                         header = T,check.names=FALSE)
homo<-read.table("homo.txt",header=T,sep='\t')
RNAseq_expr<-merge(homo,LIHC_gene2,by.x="GeneID",by.y="GeneID",all=FALSE)
RNAseq_expr2<-RNAseq_expr[,-c(1,3,4)]
RNAseq_expr3<-RNAseq_expr[,-c(1,2,3,4)]
col_RNAseq_expr<-data.frame(RNAseq_expr$Symbol)
colnames(col_RNAseq_expr)<-'gene'

tt<-data.frame(colnames(RNAseq_expr2)[-1])
colnames(tt)<-"sample"
sample_col<-separate(tt,col = "sample",into = c("v1","v2","v3","v4","v5","v6","v7"),remove = F,fill="right")

tumor_rna<-RNAseq_expr3[,sample_col$v4=="01A"]
dim(tumor_rna)
la<-data.frame(rep(0,ncol(tumor_rna)))
colnames(la)<-"if_health"

adj_tumor_rna<-RNAseq_expr3[,sample_col$v4=="11A"]
dim(adj_tumor_rna)
la2<-data.frame(rep(1,ncol(adj_tumor_rna)))
colnames(la2)<-"if_health"

label_6<-rbind(la,la2)

RNAseq_nolabel<-cbind(tumor_rna,adj_tumor_rna)
my_RNAseq<-cbind(col_RNAseq_expr,RNAseq_nolabel)



##t-test top 33 genes  &  clu3 33 genes
#gpf<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t")
gpf<-read.table("p_limma_order.txt",header = TRUE,sep="\t")

data<-data.frame(gpf$genename[1:33])
colnames(data)<-"gene_select"

clu3gene<-read.table("Genes_SelectbyRFEof_expr_label_cluster 3 .txt",header = T,sep="\t")

##all merge
same_col1<-merge(genename_pvalue_fdr,col_GSE25097_expr,by.x="genename",by.y="gene",all=FALSE)
same_col2<-merge(same_col1,col_GSE14520_expr,by.x="genename",by.y="gene",all=FALSE)
same_col3<-merge(same_col2,col_GSE22058_expr,by.x="genename",by.y="gene",all=FALSE)
same_col4<-merge(same_col3,col_GSE63898_expr,by.x="genename",by.y="gene",all=FALSE)
same_col5<-merge(same_col4,col_GSE64041_expr,by.x="genename",by.y="gene",all=FALSE)
same_col6<-merge(same_col5,col_GSE45436_expr,by.x="genename",by.y="gene",all=FALSE)
same_col7<-merge(same_col6,col_RNAseq_expr,by.x="genename",by.y="gene",all=FALSE)


diffGene1<-data.frame(setdiff(same_col7$genename,clu3gene$gene_select),check.names = F)
colnames(diffGene1)<-'genename'
diffGene2<-data.frame(setdiff(diffGene1$genename,data$gene_select),check.names = F)
colnames(diffGene2)<-'gene_select'

#7个数据集的基因表达数据
diffgene_GSE25097_expr<-merge(diffGene2,GSE25097_expr2,by.x="gene_select",by.y="gene",all=FALSE)
diffgene_GSE14520_expr<-merge(diffGene2,GSE14520_expr,by.x="gene_select",by.y="Symbol",all=FALSE)
diffgene_GSE22058_expr<-merge(diffGene2,GSE22058_expr,by.x="gene_select",by.y="Symbol",all=FALSE)
diffgene_GSE63898_expr<-merge(diffGene2,GSE63898_expr,by.x="gene_select",by.y="Symbol",all=FALSE)
diffgene_GSE64041_expr<-merge(diffGene2,GSE64041_expr,by.x="gene_select",by.y="symbol",all=FALSE)
diffgene_GSE45436_expr<-merge(diffGene2,GSE45436_expr,by.x="gene_select",by.y="Symbol",all=FALSE)
diffgene_RNAseq_expr<-merge(diffGene2,my_RNAseq,by.x="gene_select",by.y="gene",all=FALSE)


# test_fitted<-data.frame()
# test_if_health<-data.frame()
# 
# test_fitted_2<-data.frame()
# test_if_health_2<-data.frame()
# 
# test_fitted_3<-data.frame()
# test_if_health_3<-data.frame()
# 
# test_fitted_4<-data.frame()
# test_if_health_4<-data.frame()
# 
# test_fitted_5<-data.frame()
# test_if_health_5<-data.frame()
# 
# test_fitted_6<-data.frame()
# test_if_health_6<-data.frame()

My_GSE14520_auc<-c()
My_GSE22058_auc<-c()
My_GSE63898_auc<-c()
My_GSE64041_auc<-c()
My_GSE45436_auc<-c()
My_RNAseq_auc<-c()

for(j in 1:50){#
  
  Random_Gene <- data.frame(diffGene2[ sample(1:nrow(diffGene2),33,replace = F),1 ])
  colnames(Random_Gene)<-"gene_select"
  
  #GSE14520
  #test_data
  feature_gene_expr<-merge(Random_Gene,diffgene_GSE14520_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr)<-feature_gene_expr[,1]
  feature_gene_expr<-feature_gene_expr[,-1]
  a1<-apply(feature_gene_expr, 2, scale)
  sa2<-data.frame(apply(a1, 2, sigmoid))
  t_feature_gene_expr<-data.frame(t(sa2),check.names = F)
  colnames(t_feature_gene_expr)<-row.names(feature_gene_expr)
  mylabel<-gse14520_GSM_label[row.names(t_feature_gene_expr),3]
  
  feature_gene_expr_label<-cbind(t_feature_gene_expr,mylabel)
  colnames(feature_gene_expr_label)[ncol(feature_gene_expr_label)]<-"if_health"
  
  ##train_data
  GSE25097_expr_FG<-GSE25097_expr[row.names(feature_gene_expr),]
  b1<-apply(GSE25097_expr_FG, 2, scale)
  sb2<-data.frame(apply(b1, 2, sigmoid))
  t_GSE25097_expr_FG<-data.frame(t(sb2),check.names = F)
  colnames(t_GSE25097_expr_FG)<-row.names(GSE25097_expr_FG)
  train_data_label<-mutate(  t_GSE25097_expr_FG,if_health=c(rep(1,243),rep(0,268))   )
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out <- tune.svm(if_health ~ .,data=train_data_label,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out)
  bestmod<-tune.out$best.model
  
  test_pre<-predict(bestmod,feature_gene_expr_label,decision.values=TRUE)
  fitted1 <- attributes(test_pre)$decision.values
  sf<-sigmoid(fitted1)
  
  # dsf<-data.frame(sf)
  # colnames(dsf)[1]<-paste('random_',j,sep='')
  # dti<-data.frame(feature_gene_expr_label$if_health)
  # colnames(dti)[1]<-paste('random_',j,sep='')
  # 
  # if(j==1){  test_fitted<-dsf 
  # test_if_health<-dti  
  # }else{  test_fitted<-cbind(test_fitted,dsf) 
  # test_if_health<-cbind(test_if_health,dti)  }
  
  fitted<-sf  
  
  pred<-prediction(as.numeric(fitted),as.numeric(feature_gene_expr_label$if_health))
  per<-performance(pred,measure="auc")
  auc<-as.numeric(per@y.values)
  print(auc)
  My_GSE14520_auc<-c(My_GSE14520_auc,auc)
  
  # roc_clu <- plot.roc(as.numeric(feature_gene_expr_label$if_health),as.numeric(fitted),
  #                     main=paste("test data_",Clu_FileNames[j],sep=''), percent=TRUE, 
  #                     col=rainbow(60)[j])#,print.auc=T
  
  #GSE22058
  #test_data
  feature_gene_expr_2<-merge(Random_Gene,diffgene_GSE22058_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr_2)<-feature_gene_expr_2[,1]
  feature_gene_expr_2<-feature_gene_expr_2[,-1]
  a1_2<-apply(feature_gene_expr_2, 2, scale)
  sa2_2<-data.frame(apply(a1_2, 2, sigmoid))
  t_feature_gene_expr_2<-data.frame(t(sa2_2),check.names = F)
  colnames(t_feature_gene_expr_2)<-row.names(feature_gene_expr_2)
  mylabel_2<-gse22058_GSM_label[row.names(t_feature_gene_expr_2),11]
  
  feature_gene_expr_label_2<-cbind(t_feature_gene_expr_2,mylabel_2)
  colnames(feature_gene_expr_label_2)[ncol(feature_gene_expr_label_2)]<-"if_health"
  
  ##train_data
  # GSE25097_expr_FG_2<-GSE25097_expr[row.names(feature_gene_expr_2),]
  # b1_2<-apply(GSE25097_expr_FG_2, 2, scale)
  # sb2_2<-data.frame(apply(b1_2, 2, sigmoid))
  # t_GSE25097_expr_FG_2<-data.frame(t(sb2_2),check.names = F)
  # colnames(t_GSE25097_expr_FG_2)<-row.names(GSE25097_expr_FG_2)
  # train_data_label_2<-mutate(  t_GSE25097_expr_FG_2,if_health=c(rep(1,243),rep(0,268))   )
  # 
  # ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  # tune.out_2 <- tune.svm(if_health ~ .,data=train_data_label_2,kernel='radial',decision.values=TRUE,
  #                        cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  # summary(tune.out_2)
  # bestmod_2<-tune.out_2$best.model
  
  test_pre_2<-predict(bestmod,feature_gene_expr_label_2,decision.values=TRUE)
  fitted2 <- attributes(test_pre_2)$decision.values
  sf_2<-sigmoid(fitted2)
  fitted_2<-sf_2 
  
  # dsf_2<-data.frame(sf_2)
  # colnames(dsf_2)[1]<-paste('random_',j,sep='')
  # dti_2<-data.frame(feature_gene_expr_label_2$if_health)
  # colnames(dti_2)[1]<-paste('random_',j,sep='')
  # 
  # if(j==1){  test_fitted_2<-dsf_2 
  # test_if_health_2<-dti_2  
  # }else{  test_fitted_2<-cbind(test_fitted_2,dsf_2) 
  # test_if_health_2<-cbind(test_if_health_2,dti_2)  }
  # 
  
  pred2<-prediction(as.numeric(fitted_2),as.numeric(feature_gene_expr_label_2$if_health))
  per2<-performance(pred2,measure="auc")
  auc2<-as.numeric(per2@y.values)
  
  print(auc2)
  My_GSE22058_auc<-c(My_GSE22058_auc,auc2)
  
  
  ###GSE63898
  #test_data
  feature_gene_expr_3<-merge(Random_Gene,diffgene_GSE63898_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr_3)<-feature_gene_expr_3[,1]
  feature_gene_expr_3<-feature_gene_expr_3[,-1]
  a1_3<-apply(feature_gene_expr_3, 2, scale)
  sa2_3<-data.frame(apply(a1_3, 2, sigmoid))
  t_feature_gene_expr_3<-data.frame(t(sa2_3),check.names = F)
  colnames(t_feature_gene_expr_3)<-row.names(feature_gene_expr_3)
  mylabel_3<-gse63898_GSM_label[row.names(t_feature_gene_expr_3),4]
  
  feature_gene_expr_label_3<-cbind(t_feature_gene_expr_3,mylabel_3)
  colnames(feature_gene_expr_label_3)[ncol(feature_gene_expr_label_3)]<-"if_health"
  
  ##train_data
  # GSE25097_expr_FG_3<-GSE25097_expr[row.names(feature_gene_expr_2),]
  # b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
  # sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
  # t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
  # colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
  # train_data_label_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )
  # 
  # ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  # tune.out_3 <- tune.svm(if_health ~ .,data=train_data_label_3,kernel='radial',decision.values=TRUE,
  #                        cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  # summary(tune.out_3)
  # bestmod_3<-tune.out_3$best.model
  
  test_pre_3<-predict(bestmod,feature_gene_expr_label_3,decision.values=TRUE)
  fitted3 <- attributes(test_pre_3)$decision.values
  sf_3<-sigmoid(fitted3)
  fitted_3<-sf_3 
  
  # dsf_3<-data.frame(sf_3)
  # colnames(dsf_3)[1]<-paste('random_',j,sep='')
  # dti_3<-data.frame(feature_gene_expr_label_3$if_health)
  # colnames(dti_3)[1]<-paste('random_',j,sep='')
  # 
  # if(j==1){  test_fitted_3<-dsf_3 
  # test_if_health_3<-dti_3  
  # }else{  test_fitted_3<-cbind(test_fitted_3,dsf_3) 
  # test_if_health_3<-cbind(test_if_health_3,dti_3)  }
  
  
  pred3<-prediction(as.numeric(fitted_3),as.numeric(feature_gene_expr_label_3$if_health))
  per3<-performance(pred3,measure="auc")
  auc3<-as.numeric(per3@y.values)
  
  print(auc3)
  My_GSE63898_auc<-c(My_GSE63898_auc,auc3)
  
  # roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_3$if_health),as.numeric(fitted_3),
  #                     percent=TRUE, col=rainbow(60)[j+40],
  #                     add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
  
  ###GSE64041
  #test_data
  feature_gene_expr_4<-merge(Random_Gene,diffgene_GSE64041_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr_4)<-feature_gene_expr_4[,1]
  feature_gene_expr_4<-feature_gene_expr_4[,-1]
  a1_4<-apply(feature_gene_expr_4, 2, scale)
  sa2_4<-data.frame(apply(a1_4, 2, sigmoid))
  t_feature_gene_expr_4<-data.frame(t(sa2_4),check.names = F)
  colnames(t_feature_gene_expr_4)<-row.names(feature_gene_expr_4)
  
  feature_gene_expr_label_4<-mutate(t_feature_gene_expr_4,if_health=c(rep(1,5),rep(c(1,0),60)))
  #colnames(feature_gene_expr_label_4)[ncol(feature_gene_expr_label_4)]<-"if_health"
  
  ##train_data
  # GSE25097_expr_FG_3<-GSE25097_expr[row.names(feature_gene_expr_2),]
  # b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
  # sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
  # t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
  # colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
  # train_data_label_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )
  # 
  # ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  # tune.out_3 <- tune.svm(if_health ~ .,data=train_data_label_3,kernel='radial',decision.values=TRUE,
  #                        cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  # summary(tune.out_3)
  # bestmod_3<-tune.out_3$best.model
  
  test_pre_4<-predict(bestmod,feature_gene_expr_label_4,decision.values=TRUE)
  fitted4 <- attributes(test_pre_4)$decision.values
  sf_4<-sigmoid(fitted4)
  fitted_4<-sf_4 
  
  # dsf_4<-data.frame(sf_4)
  # colnames(dsf_4)[1]<-paste('random_',j,sep='')
  # dti_4<-data.frame(feature_gene_expr_label_4$if_health)
  # colnames(dti_4)[1]<-paste('random_',j,sep='')
  # 
  # if(j==1){  test_fitted_4<-dsf_4 
  # test_if_health_4<-dti_4  
  # }else{  test_fitted_4<-cbind(test_fitted_4,dsf_4) 
  # test_if_health_4<-cbind(test_if_health_4,dti_4)  }
  
  
  pred4<-prediction(as.numeric(fitted_4),as.numeric(feature_gene_expr_label_4$if_health))
  per4<-performance(pred4,measure="auc")
  auc4<-as.numeric(per4@y.values)
  
  print(auc4)
  My_GSE64041_auc<-c(My_GSE64041_auc,auc4)
  
  # roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_3$if_health),as.numeric(fitted_3),
  #                     percent=TRUE, col=rainbow(60)[j+40],
  #                     add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
  
  ###GSE45436
  #test_data
  feature_gene_expr_5<-merge(Random_Gene,diffgene_GSE45436_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr_5)<-feature_gene_expr_5[,1]
  feature_gene_expr_5<-feature_gene_expr_5[,-1]
  a1_5<-apply(feature_gene_expr_5, 2, scale)
  sa2_5<-data.frame(apply(a1_5, 2, sigmoid))
  t_feature_gene_expr_5<-data.frame(t(sa2_5),check.names = F)
  colnames(t_feature_gene_expr_5)<-row.names(feature_gene_expr_5)
  #mylabel_5<-gse45436_GSM_label[row.names(t_feature_gene_expr_5),4]
  
  feature_gene_expr_label_5<-cbind(t_feature_gene_expr_5,GSE45436_mylabel)
  colnames(feature_gene_expr_label_5)[ncol(feature_gene_expr_label_5)]<-"if_health"
  
  ##train_data
  # GSE25097_expr_FG_3<-GSE25097_expr[row.names(feature_gene_expr_2),]
  # b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
  # sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
  # t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
  # colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
  # train_data_label_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )
  # 
  # ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  # tune.out_3 <- tune.svm(if_health ~ .,data=train_data_label_3,kernel='radial',decision.values=TRUE,
  #                        cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  # summary(tune.out_3)
  # bestmod_3<-tune.out_3$best.model
  
  test_pre_5<-predict(bestmod,feature_gene_expr_label_5,decision.values=TRUE)
  fitted5 <- attributes(test_pre_5)$decision.values
  sf_5<-sigmoid(fitted5)
  fitted_5<-sf_5 
  
  # dsf_5<-data.frame(sf_5)
  # colnames(dsf_5)[1]<-paste('random_',j,sep='')
  # dti_5<-data.frame(feature_gene_expr_label_5$if_health)
  # colnames(dti_5)[1]<-paste('random_',j,sep='')
  # 
  # if(j==1){  test_fitted_5<-dsf_5 
  # test_if_health_5<-dti_5  
  # }else{  test_fitted_5<-cbind(test_fitted_5,dsf_5) 
  # test_if_health_5<-cbind(test_if_health_5,dti_5)  }
  # 
  
  pred5<-prediction(as.numeric(fitted_5),as.numeric(feature_gene_expr_label_5$if_health))
  per5<-performance(pred5,measure="auc")
  auc5<-as.numeric(per5@y.values)
  
  print(auc5)
  My_GSE45436_auc<-c(My_GSE45436_auc,auc5)
  
  # roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_3$if_health),as.numeric(fitted_3),
  #                     percent=TRUE, col=rainbow(60)[j+40],
  #                     add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
  
  
  ###RNA_seq
  #test_data
  feature_gene_expr_6<-merge(Random_Gene,diffgene_RNAseq_expr,by.x="gene_select",by.y="gene_select",all=FALSE)
  
  row.names(feature_gene_expr_6)<-feature_gene_expr_6[,1]
  feature_gene_expr_6<-feature_gene_expr_6[,-1]
  a1_6<-apply(feature_gene_expr_6, 2, scale)
  sa2_6<-data.frame(apply(a1_6, 2, sigmoid))
  t_feature_gene_expr_6<-data.frame(t(sa2_6),check.names = F)
  colnames(t_feature_gene_expr_6)<-row.names(feature_gene_expr_6)
  
  feature_gene_expr_label_6<-cbind(t_feature_gene_expr_6,label_6)
  #colnames(feature_gene_expr_label_6)[ncol(feature_gene_expr_label_6)]<-"if_health"
  
  ##train_data
  # GSE25097_expr_FG_3<-GSE25097_expr[row.names(feature_gene_expr_2),]
  # b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
  # sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
  # t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
  # colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
  # train_data_label_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )
  # 
  # ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  # tune.out_3 <- tune.svm(if_health ~ .,data=train_data_label_3,kernel='radial',decision.values=TRUE,
  #                        cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  # summary(tune.out_3)
  # bestmod_3<-tune.out_3$best.model
  
  test_pre_6<-predict(bestmod,feature_gene_expr_label_6,decision.values=TRUE)
  fitted6 <- attributes(test_pre_6)$decision.values
  sf_6<-sigmoid(fitted6)
  fitted_6<-sf_6 
  
  # dsf_6<-data.frame(sf_6)
  # colnames(dsf_6)[1]<-paste("Random_",j,sep='')
  # dti_6<-data.frame(feature_gene_expr_label_6$if_health)
  # colnames(dti_6)[1]<-paste("Random_",j,sep='')
  # 
  # if(j==1){  test_fitted_6<-dsf_6 
  # test_if_health_6<-dti_6  
  # }else{  test_fitted_6<-cbind(test_fitted_6,dsf_6) 
  # test_if_health_6<-cbind(test_if_health_6,dti_6)  }
  
  
  pred6<-prediction(as.numeric(fitted_6),as.numeric(feature_gene_expr_label_6$if_health))
  per6<-performance(pred6,measure="auc")
  auc6<-as.numeric(per6@y.values)
  
  print(auc6)
  My_RNAseq_auc<-c(My_RNAseq_auc,auc6)
  
  # roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_3$if_health),as.numeric(fitted_3),
  #                     percent=TRUE, col=rainbow(60)[j+40],
  #                     add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
  
}
mean(My_GSE14520_auc)
mean(My_GSE22058_auc)
mean(My_GSE63898_auc)
mean(My_GSE64041_auc)
mean(My_GSE45436_auc)
mean(My_RNAseq_auc)

mean_auc<-c(mean(My_GSE14520_auc),mean(My_GSE22058_auc),mean(My_GSE63898_auc),
            mean(My_GSE64041_auc),mean(My_GSE45436_auc),mean(My_RNAseq_auc))
#write.table(mean_auc,"mean_auc_random33.txt",quote=FALSE,sep="\t")
write.table(mean_auc,"limma_mean_auc_random33.txt",quote=FALSE,sep="\t")



# all_auc<-c(My_GSE14520_auc,My_GSE22058_auc,My_GSE63898_auc,My_GSE64041_auc,My_GSE45436_auc,My_RNAseq_auc)
# dim(all_auc)<-c(30,6)
# all_auc<-data.frame(all_auc)
# colnames(all_auc)<-c('My_GSE14520_auc','My_GSE22058_auc','My_GSE63898_auc','My_GSE64041_auc','My_GSE45436_auc','My_RNAseq_auc')
# row.names(all_auc)<-Clu_FileNames
# write.table(all_auc,"all_auc.txt",quote=F,sep="\t")
# 
# test_fitted<-test_fitted[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health<-test_if_health[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# test_fitted_2<-test_fitted_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health_2<-test_if_health_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# test_fitted_3<-test_fitted_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health_3<-test_if_health_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# test_fitted<-test_fitted[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health<-test_if_health[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# test_fitted_2<-test_fitted_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health_2<-test_if_health_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# test_fitted_3<-test_fitted_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# test_if_health_3<-test_if_health_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
# 
# 
# 
# write.table(test_fitted,"test_fitted.txt",quote=F,sep="\t")
# write.table(test_if_health,"test_if_health.txt",quote=F,sep="\t")
# write.table(test_fitted_2,"test_fitted_2.txt",quote=F,sep="\t")
# write.table(test_if_health_2,"test_if_health_2.txt",quote=F,sep="\t")
# write.table(test_fitted_3,"test_fitted_3.txt",quote=F,sep="\t")
# write.table(test_if_health_3,"test_if_health_3.txt",quote=F,sep="\t")
# write.table(test_fitted,"test_fitted.txt",quote=F,sep="\t")
# write.table(test_if_health,"test_if_health.txt",quote=F,sep="\t")
# write.table(test_fitted_2,"test_fitted_2.txt",quote=F,sep="\t")
# write.table(test_if_health_2,"test_if_health_2.txt",quote=F,sep="\t")
# write.table(test_fitted_3,"test_fitted_3.txt",quote=F,sep="\t")
# write.table(test_if_health_3,"test_if_health_3.txt",quote=F,sep="\t")
# 
# 
# 
