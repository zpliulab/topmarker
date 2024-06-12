setwd('D:\\A_My_Data\\HCC\\1_data\\validation_2')
library(dplyr)
library(e1071)
library(ROCR)
library(pROC)
library(ggplot2)
library(caret)
library(tidyr)
library(tidyverse)


expset3<-read.table("GSE64041_expr.txt",header=T,sep="\t")
expset4<-expset3

expset<-read.table("GSE45436_expr.txt",header=T,sep="\t")
gse45436_GSM_label<-read.table("gse45436_GSM_label.txt",header=T,sep="\t")
row.names(gse45436_GSM_label)<-gse45436_GSM_label[,1]

LIHC_gene2<-read.table("LIHC_2_gene.normalized_RNAseq__tissueTypeAll__20180815160255.txt",sep='\t',header=T)
homo<-read.table(file="homo.txt",sep="\t",header=T)


GSE25097_expr<-read.table(file="GSE25097_expr.txt",sep="\t",header=T)

Clu_Path <- "D:\\A_My_Data\\HCC\\1_data\\validation_2\\GenesSelect_byRFE"
Clu_FileNames <- dir(Clu_Path)

test_fitted<-data.frame()
test_if_health<-data.frame()

test_fitted_2<-data.frame()
test_if_health_2<-data.frame()

test_fitted_3<-data.frame()
test_if_health_3<-data.frame()

My_GSE64041_auc<-c()
My_GSE45436_auc<-c()
My_RNAseq_auc<-c()
for(j in 1:(length(Clu_FileNames))){#
  Feature_Gene<-read.table(file=paste("D:\\A_My_Data\\HCC\\1_data\\validation_2\\GenesSelect_byRFE\\",Clu_FileNames[j],sep=""),header = T)
  
  #GSE64041
  #test_data
  feature_gene_expr<-merge(Feature_Gene,expset4,by.x="gene_select",by.y = "symbol",all=FALSE)
  #setdiff(Feature_Gene$gene_select,feature_gene_expr$gene_select)#  "IGHG1" "DDOST"没有
  
  row.names(feature_gene_expr)<-feature_gene_expr[,1]
  feature_gene_expr<-feature_gene_expr[,-1]
  a1<-apply(feature_gene_expr, 2, scale)
  sa2<-data.frame(apply(a1, 2, sigmoid))
  t_feature_gene_expr<-data.frame(t(sa2),check.names = F)
  colnames(t_feature_gene_expr)<-row.names(feature_gene_expr)
  feature_gene_expr_label<-mutate(t_feature_gene_expr,if_health=c(rep(1,5),rep(c(1,0),60)))
  
  GSE25097_expr_FG<-GSE25097_expr[row.names(feature_gene_expr),]
  b1<-apply(GSE25097_expr_FG, 2, scale)
  sb2<-data.frame(apply(b1, 2, sigmoid))
  t_GSE25097_expr_FG<-data.frame(t(sb2),check.names = F)
  colnames(t_GSE25097_expr_FG)<-row.names(GSE25097_expr_FG)
  train_data_label<-mutate(t_GSE25097_expr_FG,if_health=c(rep(1,243),rep(0,268))   )
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out <- tune.svm(if_health ~ .,data=train_data_label,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out)
  bestmod<-tune.out$best.model
  
  test_pre<-predict(bestmod,feature_gene_expr_label,decision.values=TRUE)
  fitted1 <- attributes(test_pre)$decision.values
  sf<-sigmoid(fitted1)
  
  dsf<-data.frame(sf)
  colnames(dsf)[1]<-paste(Clu_FileNames[j])
  dti<-data.frame(feature_gene_expr_label$if_health)
  colnames(dti)[1]<-paste(Clu_FileNames[j])
  
  if(j==1){  test_fitted<-dsf 
             test_if_health<-dti  
  }else{  test_fitted<-cbind(test_fitted,dsf) 
         test_if_health<-cbind(test_if_health,dti)  }
  
  fitted<-sf  
  
  pred2<-prediction(as.numeric(fitted),as.numeric(feature_gene_expr_label$if_health))
  per2<-performance(pred2,measure="auc")
  auc2<-as.numeric(per2@y.values)
  print(Clu_FileNames[j])
  print(auc2)
  My_GSE64041_auc<-c(My_GSE64041_auc,auc2)
  
  roc_clu <- plot.roc(as.numeric(feature_gene_expr_label$if_health),as.numeric(fitted),
                      main=paste("test data_",Clu_FileNames[j],sep=''), percent=TRUE, 
                      col=rainbow(60)[j])#,print.auc=T
  
  #GSE45436
  #test_data
  feature_gene_expr_2<-merge(Feature_Gene,expset,by.x="gene_select",by.y = "Symbol",all=FALSE)
  #setdiff(Feature_Gene$gene_select,feature_gene_expr$gene_select)#"SIK1"  "IGHG1" "DDOST"没有
  
  row.names(feature_gene_expr_2)<-feature_gene_expr_2[,1]
  feature_gene_expr_2<-feature_gene_expr_2[,-1]
  a1_2<-apply(feature_gene_expr_2, 2, scale)
  sa2_2<-data.frame(apply(a1_2, 2, sigmoid))
  t_feature_gene_expr_2<-data.frame(t(sa2_2),check.names = F)
  colnames(t_feature_gene_expr_2)<-row.names(feature_gene_expr_2)
  
  mylabel_2<-gse45436_GSM_label[row.names(t_feature_gene_expr_2),11]
  feature_gene_expr_label_2<-cbind(t_feature_gene_expr_2,mylabel_2)
  colnames(feature_gene_expr_label_2)[ncol(feature_gene_expr_label_2)]<-"if_health"
  mylabelf_2<-data.frame(mylabel_2)
  colnames(mylabelf_2)<-"if_health"
  #write.table(mylabelf,"GSE45436_mylabel.txt",quote=F,sep='\t')
  
  
  GSE25097_expr_FG_2<-GSE25097_expr[row.names(feature_gene_expr_2),]
  b1_2<-apply(GSE25097_expr_FG_2, 2, scale)
  sb2_2<-data.frame(apply(b1_2, 2, sigmoid))
  t_GSE25097_expr_FG_2<-data.frame(t(sb2_2),check.names = F)
  colnames(t_GSE25097_expr_FG_2)<-row.names(GSE25097_expr_FG_2)
  train_data_label_2<-mutate(t_GSE25097_expr_FG_2,if_health=c(rep(1,243),rep(0,268))   )
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out_2 <- tune.svm(if_health ~ .,data=train_data_label_2,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out_2)
  bestmod_2<-tune.out_2$best.model
  
  test_pre_2<-predict(bestmod_2,feature_gene_expr_label_2,decision.values=TRUE)
  fitted1_2 <- attributes(test_pre_2)$decision.values
  sf_2<-sigmoid(fitted1_2)
  fitted_2<-sf_2  
  
  test_pre_2<-predict(bestmod_2,feature_gene_expr_label_2,decision.values=TRUE)
  fitted2_2 <- attributes(test_pre_2)$decision.values
  sf_2<-sigmoid(fitted2_2)
  fitted_2<-sf_2 
  
  dsf_2<-data.frame(sf_2)
  colnames(dsf_2)[1]<-paste(Clu_FileNames[j])
  dti_2<-data.frame(feature_gene_expr_label_2$if_health)
  colnames(dti_2)[1]<-paste(Clu_FileNames[j])
  
  if(j==1){  test_fitted_2<-dsf_2 
  test_if_health_2<-dti_2  
  }else{  test_fitted_2<-cbind(test_fitted_2,dsf_2) 
  test_if_health_2<-cbind(test_if_health_2,dti_2)  }
  
  pred3<-prediction(as.numeric(fitted_2),as.numeric(feature_gene_expr_label_2$if_health))
  per3<-performance(pred3,measure="auc")
  auc3<-as.numeric(per3@y.values)
  print(Clu_FileNames[j])
  print(auc3)
  My_GSE45436_auc<-c(My_GSE45436_auc,auc3)
  
  roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_2$if_health),as.numeric(fitted_2),
                       percent=TRUE, col=rainbow(60)[j+20],
                      add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
  
  ###RNA_Seq
  #test_data
 
  Feature_Gene_id<-merge(Feature_Gene,homo,by.x="gene_select",by.y="Symbol",all=FALSE)
  #write.table(Feature_Gene_id,"feature_gene_new_id.txt",sep='\t',quote=F)
  RNAseq_33gene<-merge(Feature_Gene_id,LIHC_gene2,by.x="GeneID",by.y="GeneID",all=FALSE)
  #write.table(RNAseq_71gene,"HCC_RNAseq_71gene.txt",sep='\t',quote=F)
  #setdiff(Feature_Gene_id$gene_select,RNAseq_71gene$gene_select)#"IGHG1"没有
  
  rna1<-data.frame(t(RNAseq_33gene[,-c(1,2,3,4)]))
  colnames(rna1)<-RNAseq_33gene$gene_select
  tt<-data.frame(row.names(rna1))
  colnames(tt)<-"sample"
  rna2<-cbind(tt,rna1)
  rna3<-separate(rna2,col = "sample",into = c("v1","v2","v3","v4","v5","v6","v7"),remove = F,fill="right")
  
  tumor_rna<-rna3[rna3$v4=="01A",-c(1,2,3,4,5,6,7,8)]
  dim(tumor_rna)
  la<-data.frame(rep(0,nrow(tumor_rna)))
  colnames(la)<-"if_health"
  t_rna<-cbind(tumor_rna,la)
  
  adj_tumor_rna<-rna3[rna3$v4=="11A",-c(1,2,3,4,5,6,7,8)]
  dim(adj_tumor_rna)
  la2<-data.frame(rep(1,nrow(adj_tumor_rna)))
  colnames(la2)<-"if_health"
  adj_t_rna<-cbind(adj_tumor_rna,la2)
  
  #setdiff(Feature_Gene_id$gene_select,LIHC_gene2$GeneName)
  rna_data<-rbind(adj_t_rna,t_rna)
  rna_data_nolabel<-rna_data[,-ncol(rna_data)]
  
  
  GSE25097_expr_FG_3<-GSE25097_expr[as.character(RNAseq_33gene$gene_select),]
  b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
  sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
  t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
  colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
  train1_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )
  
  
  trna_data_nolabel<-t(rna_data_nolabel)
  a1_3<-apply(trna_data_nolabel, 2, scale)
  sa2_3<-data.frame(apply(a1_3, 2, sigmoid))
  scale_rna_data_nolabel<-data.frame(t(sa2_3))
  colnames(scale_rna_data_nolabel)<-colnames(rna_data_nolabel)
  
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out_3 <- tune.svm(if_health ~ .,data=train1_3,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out_3)
  bestmod_3<-tune.out_3$best.model
  
  fitted_3<-c()
  for (i in 1:nrow(scale_rna_data_nolabel)) {
    tt_3<-scale_rna_data_nolabel[i,]
    test1_3<-data.frame(cbind(tt_3,rna_data$if_health[i]))
    colnames(test1_3)[dim(test1_3)[2]]<-"if_health"
    
    test_pre_3<-predict(bestmod_3,test1_3,decision.values=TRUE)
    fitted1_3 <- attributes(test_pre_3)$decision.values
    sf_3<-sigmoid(fitted1_3)
    fitted_3<-c(fitted_3,sf_3)
  }
  
  dsf_3<-data.frame(fitted_3)
  colnames(dsf_3)[1]<-paste(Clu_FileNames[j])
  dti_3<-data.frame(rna_data$if_health)
  colnames(dti_3)[1]<-paste(Clu_FileNames[j])
  
  if(j==1){  test_fitted_3<-dsf_3 
  test_if_health_3<-dti_3  
  }else{  test_fitted_3<-cbind(test_fitted_3,dsf_3) 
  test_if_health_3<-cbind(test_if_health_3,dti_3)  }
  
  pred4<-prediction(as.numeric(fitted_3),as.numeric(rna_data$if_health))
  per4<-performance(pred4,measure="auc")
  auc4<-as.numeric(per4@y.values)
  print(Clu_FileNames[j])
  print(auc4)
  My_RNAseq_auc<-c(My_RNAseq_auc,auc4)
  
  roc_clu <- plot.roc(as.numeric(rna_data$if_health),as.numeric(fitted_3),
                      percent=TRUE, col=rainbow(60)[j+40],
                      add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,
}


all_auc<-c(My_GSE64041_auc,My_GSE45436_auc,My_RNAseq_auc)
dim(all_auc)<-c(20,3)
all_auc<-data.frame(all_auc)
colnames(all_auc)<-c('My_GSE64041_auc','My_GSE45436_auc','My_RNAseq_auc')
row.names(all_auc)<-Clu_FileNames
write.table(all_auc,"all_auc_v2.txt",quote=F,sep="\t")

test_fitted<-test_fitted[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health<-test_if_health[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]

test_fitted_2<-test_fitted_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health_2<-test_if_health_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]

test_fitted_3<-test_fitted_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health_3<-test_if_health_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]




write.table(test_fitted,"v2_test_fitted.txt",quote=F,sep="\t")
write.table(test_if_health,"v2_test_if_health.txt",quote=F,sep="\t")
write.table(test_fitted_2,"v2_test_fitted_2.txt",quote=F,sep="\t")
write.table(test_if_health_2,"v2_test_if_health_2.txt",quote=F,sep="\t")
write.table(test_fitted_3,"v2_test_fitted_3.txt",quote=F,sep="\t")
write.table(test_if_health_3,"v2_test_if_health_3.txt",quote=F,sep="\t")


