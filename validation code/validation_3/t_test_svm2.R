setwd('D:\\A_My_Data\\HCC\\1_data\\validation_3')
library(dplyr)
library(e1071)
library(ROCR)
library(pROC)
library(ggplot2)
library(caret)

GSE14520_expr<-read.table("GSE14520_expr.txt",header=T,sep="\t")
gse14520_GSM_label<-read.table("gse14520_GSM_label.txt",header=T,sep="\t")
row.names(gse14520_GSM_label)<-gse14520_GSM_label[,1]

GSE22058_expr<-read.table("GSE22058_expr.txt",header=T,sep="\t")
gse22058_GSM_label<-read.table("gse22058_GSM_label.txt",header=T,sep="\t")
row.names(gse22058_GSM_label)<-gse22058_GSM_label[,1]

GSE63898_expr<-read.table("GSE63898_expr.txt",header=T,sep="\t")
gse63898_GSM_label<-read.table("gse63898_GSM_label.txt",header=T,sep="\t")
row.names(gse63898_GSM_label)<-gse63898_GSM_label[,1]

#gpf<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t")
gpf<-read.table("p_limma_order.txt",header = TRUE,sep="\t")

data<-data.frame(gpf$genename[1:33])
colnames(data)<-"gene_select"
Feature_Gene<-data

GSE25097_expr<-read.table(file="GSE25097_expr.txt",sep="\t",header=T)
# 
# Clu_Path <- "D:\\A_My_Data\\HCC\\1_data\\data_svm\\GenesSelect_byRFE"
# Clu_FileNames <- dir(Clu_Path)
# 
# test_fitted<-data.frame()
# test_if_health<-data.frame()
# 
# test_fitted_2<-data.frame()
# test_if_health_2<-data.frame()
# 
# test_fitted_3<-data.frame()
# test_if_health_3<-data.frame()
# 
# My_GSE14520_auc<-c()
# My_GSE22058_auc<-c()
# My_GSE63898_auc<-c()
# for(j in 1:(length(Clu_FileNames))){#
  

j=1
#GSE14520
#test_data
feature_gene_expr<-merge(Feature_Gene,GSE14520_expr,by.x="gene_select",by.y = "Symbol",all=FALSE)
#setdiff(Feature_Gene$gene,feature_gene_expr$gene)#"SIK1"  "IGHG1" "DDOST"没有

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
# colnames(dsf)[1]<-paste(Clu_FileNames[j])
# dti<-data.frame(feature_gene_expr_label$if_health)
# colnames(dti)[1]<-paste(Clu_FileNames[j])
# 
# if(j==1){  test_fitted<-dsf 
# test_if_health<-dti  
# }else{  test_fitted<-cbind(test_fitted,dsf) 
# test_if_health<-cbind(test_if_health,dti)  }

fitted<-sf  

pred2<-prediction(as.numeric(fitted),as.numeric(feature_gene_expr_label$if_health))
per2<-performance(pred2,measure="auc")
auc2<-as.numeric(per2@y.values)
#print(Clu_FileNames[j])
print(auc2)
#My_GSE14520_auc<-c(My_GSE14520_auc,auc2)

roc_clu <- plot.roc(as.numeric(feature_gene_expr_label$if_health),as.numeric(fitted),
                    main=paste("test data_",Clu_FileNames[j],sep=''), percent=TRUE, 
                    col=rainbow(60)[j])#,print.auc=T

#GSE22058
#test_data
feature_gene_expr_2<-merge(Feature_Gene,GSE22058_expr,by.x="gene_select",by.y = "Symbol",all=FALSE)
#setdiff(Feature_Gene$gene,feature_gene_expr_2$gene)#"SIK1"  "IGHG1" "DDOST"没有

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
GSE25097_expr_FG_2<-GSE25097_expr[row.names(feature_gene_expr_2),]
b1_2<-apply(GSE25097_expr_FG_2, 2, scale)
sb2_2<-data.frame(apply(b1_2, 2, sigmoid))
t_GSE25097_expr_FG_2<-data.frame(t(sb2_2),check.names = F)
colnames(t_GSE25097_expr_FG_2)<-row.names(GSE25097_expr_FG_2)
train_data_label_2<-mutate(  t_GSE25097_expr_FG_2,if_health=c(rep(1,243),rep(0,268))   )

###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
tune.out_2 <- tune.svm(if_health ~ .,data=train_data_label_2,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out_2)
bestmod_2<-tune.out_2$best.model

test_pre_2<-predict(bestmod_2,feature_gene_expr_label_2,decision.values=TRUE)
fitted2 <- attributes(test_pre_2)$decision.values
sf_2<-sigmoid(fitted2)
fitted_2<-sf_2 

# dsf_2<-data.frame(sf_2)
# colnames(dsf_2)[1]<-paste(Clu_FileNames[j])
# dti_2<-data.frame(feature_gene_expr_label_2$if_health)
# colnames(dti_2)[1]<-paste(Clu_FileNames[j])
# 
# if(j==1){  test_fitted_2<-dsf_2 
# test_if_health_2<-dti_2  
# }else{  test_fitted_2<-cbind(test_fitted_2,dsf_2) 
# test_if_health_2<-cbind(test_if_health_2,dti_2)  }
# 

pred3<-prediction(as.numeric(fitted_2),as.numeric(feature_gene_expr_label_2$if_health))
per3<-performance(pred3,measure="auc")
auc3<-as.numeric(per3@y.values)
print(Clu_FileNames[j])
print(auc3)
#My_GSE22058_auc<-c(My_GSE22058_auc,auc3)

roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_2$if_health),as.numeric(fitted_2),
                    percent=TRUE, col=rainbow(60)[j+20],
                    add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,

###GSE63898
#test_data
feature_gene_expr_3<-merge(Feature_Gene,GSE63898_expr,by.x="gene_select",by.y = "Symbol",all=FALSE)
#setdiff(Feature_Gene$gene,feature_gene_expr_2$gene)#"SIK1"  "IGHG1" "DDOST"没有

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
GSE25097_expr_FG_3<-GSE25097_expr[row.names(feature_gene_expr_2),]
b1_3<-apply(GSE25097_expr_FG_3, 2, scale)
sb2_3<-data.frame(apply(b1_3, 2, sigmoid))
t_GSE25097_expr_FG_3<-data.frame(t(sb2_3),check.names = F)
colnames(t_GSE25097_expr_FG_3)<-row.names(GSE25097_expr_FG_3)
train_data_label_3<-mutate(  t_GSE25097_expr_FG_3,if_health=c(rep(1,243),rep(0,268))   )

###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
tune.out_3 <- tune.svm(if_health ~ .,data=train_data_label_3,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out_3)
bestmod_3<-tune.out_3$best.model

test_pre_3<-predict(bestmod_3,feature_gene_expr_label_3,decision.values=TRUE)
fitted3 <- attributes(test_pre_3)$decision.values
sf_3<-sigmoid(fitted3)
fitted_3<-sf_3 

# dsf_3<-data.frame(sf_3)
# colnames(dsf_3)[1]<-paste(Clu_FileNames[j])
# dti_3<-data.frame(feature_gene_expr_label_3$if_health)
# colnames(dti_3)[1]<-paste(Clu_FileNames[j])
# 
# if(j==1){  test_fitted_3<-dsf_3 
# test_if_health_3<-dti_3  
# }else{  test_fitted_3<-cbind(test_fitted_3,dsf_3) 
# test_if_health_3<-cbind(test_if_health_3,dti_3)  }


pred4<-prediction(as.numeric(fitted_3),as.numeric(feature_gene_expr_label_3$if_health))
per4<-performance(pred4,measure="auc")
auc4<-as.numeric(per4@y.values)
#print(Clu_FileNames[j])
print(auc4)
# My_GSE63898_auc<-c(My_GSE63898_auc,auc4)

roc_clu <- plot.roc(as.numeric(feature_gene_expr_label_3$if_health),as.numeric(fitted_3),
                    percent=TRUE, col=rainbow(60)[j+40],
                    add=TRUE)#main=paste("test data_GSE22058 Data_",Clu_FileNames[j],sep=''),print.auc=T,

nn<-c(auc2,auc3,auc4)


all_auc<-c(My_GSE14520_auc,My_GSE22058_auc,My_GSE63898_auc)
dim(all_auc)<-c(20,3)
all_auc<-data.frame(all_auc)
colnames(all_auc)<-c('My_GSE14520_auc','My_GSE22058_auc','My_GSE63898_auc')
row.names(all_auc)<-Clu_FileNames
write.table(all_auc,"all_auc.txt",quote=F,sep="\t")

test_fitted<-test_fitted[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health<-test_if_health[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]

test_fitted_2<-test_fitted_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health_2<-test_if_health_2[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]

test_fitted_3<-test_fitted_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
test_if_health_3<-test_if_health_3[,c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]




write.table(test_fitted,"test_fitted.txt",quote=F,sep="\t")
write.table(test_if_health,"test_if_health.txt",quote=F,sep="\t")
write.table(test_fitted_2,"test_fitted_2.txt",quote=F,sep="\t")
write.table(test_if_health_2,"test_if_health_2.txt",quote=F,sep="\t")
write.table(test_fitted_3,"test_fitted_3.txt",quote=F,sep="\t")
write.table(test_if_health_3,"test_if_health_3.txt",quote=F,sep="\t")



#####
clu2<-read.table(file=paste("D:\\A_My_Data\\HCC\\1_data\\data_svm\\GenesSelect_byRFE\\",Clu_FileNames[12],sep=""),header = T)
clu3<-read.table(file=paste("D:\\A_My_Data\\HCC\\1_data\\data_svm\\GenesSelect_byRFE\\",Clu_FileNames[14],sep=""),header = T)
clu10<-read.table(file=paste("D:\\A_My_Data\\HCC\\1_data\\data_svm\\GenesSelect_byRFE\\",Clu_FileNames[2],sep=""),header = T)



feature_gene_new<-rbind(clu2,clu3,clu10)
write.table(feature_gene_new,"feature_gene_new.txt",quote=F,sep="\t")

######
#k-fold svm
#train_data3<-train_data2[,as.character(colnames(train_data))]
expr<-feature_gene_expr_label
flds <- createFolds(1:nrow(expr), k = 10, list = TRUE, returnTrain = FALSE)
all_seq<-c()
ten_fitted<-c()
for (i in 1:10){
  test_data_2<-expr[flds[[i]],]
  all_seq<-c(all_seq,flds[[i]])
  train_data_2<-expr[-flds[[i]],]
  
  svmfit = svm(if_health ~ ., train_data_2, type='C-classification', kernel='radial',
               decision.values=TRUE)
  pred = predict(svmfit, test_data_2, decision.values=TRUE)
  fitted <- attributes(pred)$decision.values
  aa<-sigmoid(fitted)
  ten_fitted<-c(ten_fitted,aa)
}

# tune.out <- tune.svm(if_health ~ .,data=train_data,kernel='radial',gamma=0.1,cost=10,
#                      cost=c(0.1,0.5,1,10,100),gamma=c(0.1,0.5,1,2,3,4,10))
# summary(tune.out)
# pred <- predict(tune.out$best.model,newx=dat[-train,])
# table(true=dat[-train,]$y,pred)

all_seq_data<-expr[all_seq,]
svmf<-data.frame(ten_fitted)
pred2<-prediction(svmf,all_seq_data$if_health)
b<-performance(pred2,measure="auc")
auc<-as.numeric(b@y.values)
print(auc)

roc_clu <- plot.roc(as.numeric(all_seq_data$if_health),as.numeric(ten_fitted),
                    main=paste("ROC For 10-fold"), percent=TRUE, col='blue',print.auc=T)










#对其他数据中的rfe筛选出的基因用svm检验分类效果。

gse_14520<-read.table("gse14520_expr.txt",sep="\t",header = T)
gse_14520_col<-data.frame(colnames(gse_14520))
colnames(gse_14520_col)<-"Name"

gse_63898<-read.table("gse63898_expr.txt",sep="\t",header = T)
gse_63898_col<-data.frame(colnames(gse_63898))
colnames(gse_63898_col)<-"Name"

gse_22058<-read.table("gse22058_expr.txt",sep="\t",header = T)
gse_22058_col<-data.frame(colnames(gse_22058))
colnames(gse_22058_col)<-"Name"

BigClu_svmData <- read.table("gse25097_expr.txt",sep="\t",header = T)

#####
Clu_Path <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\ClusterData"
Clu_FileNames <- dir(Clu_Path)

My_AUC_clu<-c()
for(j in 1:(length(Clu_FileNames))){
  Feature_Gene <- read.table(file = paste(Clu_FileNames[j],sep = "\t"),header = T)
  
  col_inter1<-merge(Feature_Gene,gse_14520_col,by.x="gene",by.y="Name",all = FALSE)
  col_inter2<-merge(col_inter1,gse_63898_col,by.x="gene",by.y="Name",all = FALSE)
  col_inter3<-merge(col_inter2,gse_22058_col,by.x="gene",by.y="Name",all = FALSE)
  
  #合并三个数据作为test_data
  gene_overlap<-as.character(col_inter3$gene)
  gse_14520_gene_overlap<-gse_14520[,gene_overlap]
  gse_63898_gene_overlap<-gse_63898[,gene_overlap]
  gse_22058_gene_overlap<-gse_22058[,gene_overlap]
  gse_data_nolabel<-rbind(gse_14520_gene_overlap,gse_63898_gene_overlap,gse_22058_gene_overlap)
  gse_label<-c(gse_14520$if_health,gse_63898$if_health,gse_22058$if_health)
  gse_labelf<-data.frame(as.factor(gse_label))
  colnames(gse_labelf)<-"if_health"
  gse_data<-cbind(gse_data_nolabel,gse_labelf)
  
  #BigClu_svmData数据处理作为train_data
  train_data_nolabel <- BigClu_svmData[,as.character(col_inter3$gene)]
  train_data <- cbind(train_data_nolabel,as.factor(BigClu_svmData$if_health))#原始数据做为训练集
  colnames(train_data)[dim(train_data)[2]] <- "if_health"
  
  tgse_data_nolabel<-t(gse_data_nolabel)
  a1<-apply(tgse_data_nolabel, 2, scale)
  sa2<-data.frame(apply(a1, 2, sigmoid))
  scale_gse_data_nolabel<-data.frame(t(sa2))
  colnames(scale_gse_data_nolabel)<-colnames(gse_data_nolabel)
 
  ttrain_data_nolabel<-t(train_data_nolabel)
  b1<-apply(ttrain_data_nolabel, 2, scale)
  sb2<-data.frame(apply(b1, 2, sigmoid))
  b2<-data.frame(t(sb2))
  train1<-data.frame(cbind(b2,train_data$if_health))
  colnames(train1)<-colnames(train_data)
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out <- tune.svm(if_health ~ .,data=train1,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out)
  bestmod<-tune.out$best.model
  
  fitted<-c()
  for (i in 1:nrow(scale_gse_data_nolabel)) {
    tt<-scale_gse_data_nolabel[i,]
    test1<-data.frame(cbind(tt,gse_data$if_health[i]))
    colnames(test1)[dim(test1)[2]]<-"if_health"
    
    test_pre<-predict(bestmod,test1,decision.values=TRUE)
    fitted1 <- attributes(test_pre)$decision.values
    fitted<-c(fitted,fitted1)
  }
  
  pred2<-prediction(as.numeric(fitted),as.numeric(gse_data$if_health))
  per2<-performance(pred2,measure="auc")
  auc2<-as.numeric(per2@y.values)
  
  My_AUC_clu<-c(My_AUC_clu,auc2)
  
  roc_clu <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted),
                      main=paste(Clu_FileNames[j]), percent=TRUE, col="blue",print.auc=T)
}  

My_AUC_cluf<-data.frame(My_AUC_clu)
write.table(My_AUC_cluf,"My_AUC_clu.txt",quote = F,sep='\t')

######
#GenesSelect_byRFE
RFE_Path <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\GenesSelect_byRFE"
RFE_FileNames <- dir(RFE_Path)

My_AUC<-c()
fitted<-data.frame()
for(k in 1:(length(RFE_FileNames))){#
  Feature_Gene <- read.table(file = paste(RFE_FileNames[k],sep = "\t"),header = T)
  col_inter1<-merge(Feature_Gene,gse_14520_col,by.x="gene_select",by.y="Name",all = FALSE)
  col_inter2<-merge(col_inter1,gse_63898_col,by.x="gene_select",by.y="Name",all = FALSE)
  col_inter3<-merge(col_inter2,gse_22058_col,by.x="gene_select",by.y="Name",all = FALSE)
  
  #合并三个数据作为test_data
  gene_overlap<-as.character(col_inter3$gene)
  gse_14520_gene_overlap<-gse_14520[,gene_overlap]
  gse_63898_gene_overlap<-gse_63898[,gene_overlap]
  gse_22058_gene_overlap<-gse_22058[,gene_overlap]
  gse_data_nolabel<-rbind(gse_14520_gene_overlap,gse_63898_gene_overlap,gse_22058_gene_overlap)
  gse_label<-c(gse_14520$if_health,gse_63898$if_health,gse_22058$if_health)
  gse_labelf<-data.frame(as.factor(gse_label))
  colnames(gse_labelf)<-"if_health"
  gse_data<-cbind(gse_data_nolabel,gse_labelf)
  
  #BigClu_svmData数据处理作为train_data
  train_data_nolabel <- BigClu_svmData[,as.character(col_inter3$gene)]
  train_data <- cbind(train_data_nolabel,as.factor(BigClu_svmData$if_health))#原始数据做为训练集
  colnames(train_data)[dim(train_data)[2]] <- "if_health"
  
  tgse_data_nolabel<-t(gse_data_nolabel)
  a1<-apply(tgse_data_nolabel, 2, scale)
  sa2<-data.frame(apply(a1, 2, sigmoid))
  scale_gse_data_nolabel<-data.frame(t(sa2))
  colnames(scale_gse_data_nolabel)<-colnames(gse_data_nolabel)
  # 
  ttrain_data_nolabel<-t(train_data_nolabel)
  b1<-apply(ttrain_data_nolabel, 2, scale)
  sb2<-data.frame(apply(b1, 2, sigmoid))
  b2<-data.frame(t(sb2))
  train1<-data.frame(cbind(b2,train_data$if_health))
  colnames(train1)<-colnames(train_data)

  #heatmap(as.matrix(gse_data_nolabel), scale = "none", col = bluered(100))
  
  ###合并一个待检测样本到训练样本中，进行归一化，再进行svm分类。
  tune.out <- tune.svm(if_health ~ .,data=train1,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))#,0.5,1,10,scale=TRUE
  summary(tune.out)
  bestmod<-tune.out$best.model
  
  
  for (i in 1:nrow(scale_gse_data_nolabel)) {
    tt<-scale_gse_data_nolabel[i,]
    test1<-data.frame(cbind(tt,gse_data$if_health[i]))
    colnames(test1)[dim(test1)[2]]<-"if_health"
    
    test_pre<-predict(bestmod,test1,decision.values=TRUE)
    fitted1 <- attributes(test_pre)$decision.values
    fitted[i,k]<-fitted1
  }
  
  pred2<-prediction(as.numeric(fitted[,k]),as.numeric(gse_data$if_health))
  b2<-performance(pred2,measure="auc")
  auc2<-as.numeric(b2@y.values)
  My_AUC<-c(My_AUC,auc2)
  
  roc_clu <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,k]),
                      main=paste(RFE_FileNames[k]), percent=TRUE, col="blue",print.auc=T)
  
} 

My_AUCf<-data.frame(My_AUC)
write.table(My_AUCf,"My_AUC.txt",quote = F,sep='\t')



roc_clu1 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,1]), col=rainbow(20)[1]) 
roc_clu10 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,2]), col=rainbow(20)[2])
roc_clu11 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,3]), col=rainbow(20)[3])
roc_clu12 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,4]), col=rainbow(20)[4])
roc_clu13 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,5]), col=rainbow(20)[5])
roc_clu14 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,6]), col='black')
roc_clu15 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,7]), col=rainbow(20)[7])
roc_clu16 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,8]), col=rainbow(20)[8])
roc_clu17 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,9]), col=rainbow(20)[9])
roc_clu18 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,10]), col=rainbow(20)[10])
roc_clu19 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,11]), col=rainbow(20)[11])
roc_clu2 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,12]), col='#470024')
roc_clu20 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,13]), col=rainbow(20)[13])
roc_clu3 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,14]), col='blue')
roc_clu4 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,15]), col=rainbow(20)[15])
roc_clu5 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,16]), col=rainbow(20)[16])
roc_clu6 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,17]), col=rainbow(20)[17])
roc_clu7 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,18]), col=rainbow(20)[18])
roc_clu8 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,19]), col=rainbow(20)[19])
roc_clu9 <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted[,20]), col=rainbow(20)[20])


all_roc<-ggroc(list(cluster1=roc_clu1,cluster2=roc_clu2,cluster3=roc_clu3,cluster4=roc_clu4,
                    cluster5=roc_clu5,cluster6=roc_clu6,cluster7=roc_clu7,cluster8=roc_clu8,
                    cluster9=roc_clu9,cluster10=roc_clu10,cluster11=roc_clu11,cluster12=roc_clu12,
                    cluster13=roc_clu13,cluster14=roc_clu14,cluster15=roc_clu15,cluster16=roc_clu16,
                    cluster17=roc_clu17,cluster18=roc_clu18,cluster19=roc_clu19,cluster20=roc_clu20),legacy.axes=TRUE)
plot(all_roc)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
all_roc + ggtitle("My ROC curve") + scale_fill_manual(values=cbbPalette) 


######
#GenesSelect_byRFE_FeatureData
RFE_FD_Path <- "C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\GenesSelect_byRFE_FeatureData"
RFE_FD_FileNames <- dir(RFE_FD_Path)


# name<-data.frame(c(paste("cluster",1:20,sep="")))
# colnames(name)<-"name"
#MY_ROC<-c()
#my_list<-list()
My_AUC<-c()
for(k in 2:(length(RFE_FD_FileNames))){#
  Feature_Gene <- read.table(file = paste(RFE_FD_FileNames[k],sep = "\t"),header = T)
  col_inter1<-merge(Feature_Gene,gse_14520_col,by.x="gene_select",by.y="Name",all = FALSE)
  col_inter2<-merge(col_inter1,gse_63898_col,by.x="gene_select",by.y="Name",all = FALSE)
  col_inter3<-merge(col_inter2,gse_22058_col,by.x="gene_select",by.y="Name",all = FALSE)
  
  #合并三个数据作为test_data
  gene_overlap<-as.character(col_inter3$gene)
  gse_14520_gene_overlap<-gse_14520[,gene_overlap]
  gse_63898_gene_overlap<-gse_63898[,gene_overlap]
  gse_22058_gene_overlap<-gse_22058[,gene_overlap]
  gse_data_nolabel<-rbind(gse_14520_gene_overlap,gse_63898_gene_overlap,gse_22058_gene_overlap)
  gse_label<-c(gse_14520$if_health,gse_63898$if_health,gse_22058$if_health)
  gse_labelf<-data.frame(as.factor(gse_label))
  colnames(gse_labelf)<-"if_health"
  gse_data<-cbind(gse_data_nolabel,gse_labelf)
  
  #BigClu_svmData数据处理作为train_data
  train_data_nolabel <- BigClu_svmData[,as.character(col_inter3$gene)]
  train_data <- cbind(train_data_nolabel,as.factor(BigClu_svmData$if_health))#原始数据做为训练集
  colnames(train_data)[dim(train_data)[2]] <- "if_health"
  
  tgse_data_nolabel<-t(gse_data_nolabel)
  a1<-apply(tgse_data_nolabel, 2, scale)
  sa2<-data.frame(apply(a1, 2, sigmoid))
  scale_gse_data_nolabel<-data.frame(t(sa2))
  colnames(scale_gse_data_nolabel)<-colnames(gse_data_nolabel)
  # 
  # b1<-apply(train_data_nolabel, 2, scale)
  # b2<-data.frame(apply(b1, 2, sigmoid))
  # train1<-data.frame(cbind(b2,train_data$if_health))
  # colnames(train1)[dim(train1)[2]]<-"if_health"
  
  ttrain_data_nolabel<-t(train_data_nolabel)
  b1<-apply(ttrain_data_nolabel, 2, scale)
  sb2<-data.frame(apply(b1, 2, sigmoid))
  b2<-data.frame(t(sb2))
  train1<-data.frame(cbind(b2,train_data$if_health))
  colnames(train1)<-colnames(train_data)
  
  #heatmap(as.matrix(gse_data_nolabel), scale = "none", col = bluered(100))

  tune.out <- tune.svm(if_health ~ .,data=train1,kernel='radial',decision.values=TRUE,
                       cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
  summary(tune.out)
  bestmod<-tune.out$best.model
  
  fitted<-c()
  for (i in 1:nrow(scale_gse_data_nolabel)) {
    tt<-scale_gse_data_nolabel[i,]
    test1<-data.frame(cbind(tt,gse_data$if_health[i]))
    colnames(test1)[dim(test1)[2]]<-"if_health"
    
    test_pre<-predict(bestmod,test1,decision.values=TRUE)
    fitted1 <- attributes(test_pre)$decision.values
    fitted<-c(fitted,fitted1)
    
  }
  
  pred2<-prediction(as.numeric(fitted),as.numeric(gse_data$if_health))
 
  b2<-performance(pred2,measure="auc")
  auc2<-as.numeric(b2@y.values)
  My_AUC<-c(My_AUC,auc2)
  
  roc_clu <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted),
                      main=paste(RFE_FD_FileNames[k]), percent=TRUE, col="blue",print.auc=T)
  
  # roc_clu <- plot.roc(as.numeric(gse_data$if_health),as.numeric(fitted), col=rainbow(20)[k])
  # 
  # my_list<-list(my_list,assign(paste("cluster",k,sep=""),roc_clu))
  # 
  
}  


# all_roc<-ggroc(my_list)
# plot(all_roc)  





#####
#random genes
gse14520_all_expr <- read.table("gse14520_all_expr.txt",sep="\t",header = T)
gse14520_all_expr_col<-data.frame(colnames(gse14520_all_expr)[-dim(gse14520_all_expr)[2]])
colnames(gse14520_all_expr_col)<-"Name"

gse63898_all_expr <- read.table("gse63898_all_expr.txt",sep="\t",header = T)
gse63898_all_expr_col<-data.frame(colnames(gse63898_all_expr)[-dim(gse63898_all_expr)[2]])
colnames(gse63898_all_expr_col)<-"Name"

gse22058_all_expr <- read.table("gse22058_all_expr.txt",sep="\t",header = T)
gse22058_all_expr_col<-data.frame(colnames(gse22058_all_expr)[-dim(gse22058_all_expr)[2]])
colnames(gse22058_all_expr_col)<-"Name"

myData_all_expr <- read.table("gse25097_all_expr.txt",sep="\t",header = T)
myData_all_expr_col<-data.frame(colnames(myData_all_expr)[-dim(myData_all_expr)[2]])
colnames(myData_all_expr_col)<-"Name"

col_overlap1<-merge(gse14520_all_expr_col,gse63898_all_expr_col,by.x="Name",by.y="Name",all = FALSE)
col_overlap2<-merge(col_overlap1,gse22058_all_expr_col,by.x="Name",by.y="Name",all = FALSE)
col_overlap3<-merge(col_overlap2,myData_all_expr_col,by.x="Name",by.y="Name",all = FALSE)

genename_pvalue_fdr<-read.table("genename_pvalue_fdr.txt",header = TRUE,sep="\t")
col_overlap4<-merge(col_overlap3,genename_pvalue_fdr,by.x="Name",by.y="genename",all = FALSE)
#aa<-merge(myData_all_expr_col,genename_pvalue_fdr,by.x="Name",by.y="genename",all = FALSE)

Random_Gene <- data.frame(col_overlap4[sample(1:1600,71,replace = F),1])
Random_Gene <- data.frame(col_overlap4[sample(1601:3200,71,replace = F),1])
Random_Gene <- data.frame(col_overlap4[sample(3201:4800,71,replace = F),1])
Random_Gene <- data.frame(col_overlap4[sample(4801:6400,71,replace = F),1])
Random_Gene <- data.frame(col_overlap4[sample(6401:8219,71,replace = F),1])
colnames(Random_Gene)<-"gene"
# Random_Gene <- data.frame(col_overlap3[sample(nrow(col_overlap3),36,replace = F),1])


#合并三个数据作为test_data
gene_overlap<-as.character(Random_Gene$gene)
gse_14520_gene_overlap<-gse14520_all_expr[,gene_overlap]
gse_63898_gene_overlap<-gse63898_all_expr[,gene_overlap]
gse_22058_gene_overlap<-gse22058_all_expr[,gene_overlap]
gse_data_nolabel<-rbind(gse_14520_gene_overlap,gse_63898_gene_overlap,gse_22058_gene_overlap)
gse_label<-c(gse14520_all_expr$if_health,gse63898_all_expr$if_health,gse22058_all_expr$if_health)
gse_labelf<-data.frame(as.factor(gse_label))
colnames(gse_labelf)<-"if_health"
gse_data<-cbind(gse_data_nolabel,gse_labelf)
# write.table(gse_data,"rfe_test_expr_label.txt",quote=F,sep='\t')

#myData_svmData数据处理作为train_data
train_data_nolabel <- myData_all_expr[,as.character(Random_Gene$gene)]
train_data <- cbind(train_data_nolabel,as.factor(myData_all_expr$if_health))#原始数据做为训练集
colnames(train_data)[dim(train_data)[2]] <- "if_health"


#heatmap(as.matrix(gse_data_nolabel), scale = "none", col = bluered(100))

tgse_data_nolabel<-t(gse_data_nolabel)
a1<-apply(tgse_data_nolabel, 2, scale)
sa2<-data.frame(apply(a1, 2, sigmoid))
scale_gse_data_nolabel<-data.frame(t(sa2))
colnames(scale_gse_data_nolabel)<-colnames(gse_data_nolabel)
# 
ttrain_data_nolabel<-t(train_data_nolabel)
b1<-apply(ttrain_data_nolabel, 2, scale)
sb2<-data.frame(apply(b1, 2, sigmoid))
b2<-data.frame(t(sb2))
train1<-data.frame(cbind(b2,train_data$if_health))
colnames(train1)<-colnames(train_data)


tune.out <- tune.svm(if_health ~ .,data=train1,kernel='radial',decision.values=TRUE,
                                         cost=c(0.1,0.5,1,5,10),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out)
bestmod<-tune.out$best.model

fitted<-c()
test1<-cbind(scale_gse_data_nolabel,gse_data$if_health)
try1<-test1[1:400,]
svmfit<- svm(if_health ~ .,data=train1,kernel='radial',decision.values=TRUE)
                     
test_pre<-predict(svmfit,try1,decision.values=TRUE)
fitted1 <- attributes(test_pre)$decision.values
ttt<-sigmoid(fitted1)
fitted<-c(fitted,ttt)




fitted<-c()
for (i in 1:nrow(scale_gse_data_nolabel)) {
  tt<-scale_gse_data_nolabel[i,]
  test1<-cbind(tt,gse_data$if_health[i])
  colnames(test1)[dim(test1)[2]]<-"if_health"
 
  test_pre<-predict(bestmod,test1,decision.values=TRUE)
  fitted1 <- attributes(test_pre)$decision.values
  ttt<-sigmoid(fitted1)
  fitted<-c(fitted,ttt)
}

pred2<-prediction(as.numeric(fitted),as.numeric(gse_data$if_health))
perf2<-performance(pred2,"tpr","fpr")
plot(perf2,col=3)
points(c(0,1),c(0,1),type="l",lty=2)

b2<-performance(pred2,measure="auc")
auc2<-as.numeric(b2@y.values)
auc2

roc_clu <- plot.roc(as.numeric(gse_data[1:400,]$if_health),as.numeric(fitted),
                    main=("ROC For gse_data"), percent=TRUE, col="blue",print.auc=T)
#table(true=gse_data$if_health,pred=fitted)

  
 





######
#gse14520数据作为test
gse_14520<-read.table("gse14520_expr.txt",sep="\t",header = T)
gse_14520_col<-data.frame(colnames(gse_14520))
colnames(gse_14520_col)<-"Name"
col_inter1<-merge(Feature_Gene,gse_14520_col,by.x="gene",by.y="Name",all = FALSE)

#gse63898数据作为test
gse_63898<-read.table("gse63898_expr.txt",sep="\t",header = T)
gse_63898_col<-data.frame(colnames(gse_63898))
colnames(gse_63898_col)<-"Name"
col_inter2<-merge(col_inter1,gse_63898_col,by.x="gene",by.y="Name",all = FALSE)

#gse22058数据作为test
gse_22058<-read.table("gse22058_expr.txt",sep="\t",header = T)
gse_22058_col<-data.frame(colnames(gse_22058))
colnames(gse_22058_col)<-"Name"
col_inter3<-merge(col_inter2,gse_22058_col,by.x="gene",by.y="Name",all = FALSE)

write.table(col_inter3,"col_inter3_overlap.txt",quote=F,sep='\t',row.names = F)

#合并三个数据作为test_data
gene_overlap<-as.character(col_inter3$gene)
gse_14520_gene_overlap<-gse_14520[,gene_overlap]
gse_63898_gene_overlap<-gse_63898[,gene_overlap]
gse_22058_gene_overlap<-gse_22058[,gene_overlap]
gse_data_nolabel<-rbind(gse_14520_gene_overlap,gse_63898_gene_overlap,gse_22058_gene_overlap)
gse_label<-c(gse_14520$if_health,gse_63898$if_health,gse_22058$if_health)
gse_labelf<-data.frame(as.factor(gse_label))
colnames(gse_labelf)<-"if_health"
gse_data<-cbind(gse_data_nolabel,gse_labelf)
write.table(gse_data,"rfe_test_expr_label.txt",quote=F,sep='\t')
#三个数据分别整理成data+label
data_14520<-gse_data[1:445,]
data_63898<-gse_data[446:841,]
data_22058<-gse_data[842:1038,]

#rfe_svmData数据处理作为train_data
rfe_svmData <- read.table("gse25097_expr.txt",sep="\t",header = T)
train_data_nolabel <- rfe_svmData[,as.character(col_inter3$gene)]
train_data <- cbind(train_data_nolabel,as.factor(rfe_svmData$if_health))#原始数据做为训练集
colnames(train_data)[dim(train_data)[2]] <- "if_health"
write.table(train_data,"rfe_train_expr_label.txt",quote=F,sep='\t')


######
aa<-gse_data[842:1038,]
tt<-gse_data[1:445,]

tune.out <- tune.svm(if_health ~ .,data=aa,kernel='radial',decision.values=TRUE,
                      cost=c(0.1,0.5,1,10,100),gamma=c(0.0001,0.001,0.01,0.1),scale=TRUE)
summary(tune.out)
#bestmod<-tune.out$best.model
svmfit <- svm(if_health ~ .,data=aa,decision.values=TRUE,gamma=0.0001,cost=100)
summary(svmfit)
table(svmfit$fitted,aa$if_health)

gse_pre<-predict(svmfit,tt,decision.values=TRUE)
fitted <- attributes(gse_pre)$decision.values
roc_clu <- plot.roc(as.numeric(tt$if_health),as.numeric(fitted),
                    main=("ROC For gse"), percent=TRUE, col="red",print.auc=T)
table(true=tt$if_health,pred=gse_pre)

pred2<-prediction(as.numeric(fitted),as.numeric(gse_data$if_health))
perf2<-performance(pred2,"tpr","fpr")
plot(perf2,col=3)
points(c(0,1),c(0,1),type="l",lty=2)

b2<-performance(pred2,measure="auc")
auc2<-as.numeric(b2@y.values)
auc2




###########
#gse14520
expr1<-data_14520
set.seed(12)
training1 <- expr1[sample(nrow(expr1),0.8*nrow(expr1),replace=F),]
set.seed(12)
testing1 <- expr1[-sample(nrow(expr1),0.8*nrow(expr1),replace=F),]

tune.out <- tune.svm(if_health ~ .,data=training1,kernel='radial',decision.values=TRUE,
                     cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out)
#bestmod<-tune.out$best.model
svmfit <- svm(if_health ~ .,data=training1,decision.values=TRUE,gamma=0.0001,cost=10)
summary(svmfit)
test_pre<-predict(svmfit,testing1,decision.values=TRUE)
fitted <- attributes(test_pre)$decision.values
roc_clu <- plot.roc(as.numeric(testing1$if_health),as.numeric(fitted),
                    main=("ROC For GSE14520"), percent=TRUE, col="blue",print.auc=T)
table(true=testing1$if_health,pred=test_pre)

##random genes
gse14520_all_expr <- read.table("gse14520_all_expr.txt",sep="\t",header = T)
sub_334<-subset(gse14520_all_expr,select=-c(col_inter3$gene))
ran_gene<-sub_334[,sample((ncol(sub_334)-1),10,replace = F)]
ran_data<-cbind(ran_gene,as.factor(gse14520_all_expr$if_health))
colnames(ran_data)[ncol(ran_data)]<-"if_health"

set.seed(12)
ran_training1 <- ran_data[sample(nrow(ran_data),0.8*nrow(ran_data),replace=F),]
set.seed(12)
ran_testing1 <- ran_data[-sample(nrow(ran_data),0.8*nrow(ran_data),replace=F),]

ran_tune.out <- tune.svm(if_health ~ .,data=ran_training1,kernel='radial',decision.values=TRUE,
                     cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
#summary(ran_tune.out)
ran_bestmod<-ran_tune.out$best.model
#ran_svmfit <- svm(if_health ~ .,data=ran_training1,decision.values=TRUE,gamma=0.0001,cost=5)
#summary(ran_svmfit)
ran_test_pre<-predict(ran_bestmod,ran_testing1,decision.values=TRUE)
ran_fitted <- attributes(ran_test_pre)$decision.values
ran_roc_clu <- plot.roc(as.numeric(ran_testing1$if_health),as.numeric(ran_fitted),
                    main=("ROC For GSE14520 random genes"), percent=TRUE, col="blue",print.auc=T)
table(true=ran_testing1$if_health,pred=ran_test_pre)


#gse63898
##############
expr2<-data_63898
set.seed(123)
training2 <- expr2[sample(nrow(expr2),0.8*nrow(expr2),replace=F),]
set.seed(123)
testing2 <- expr2[-sample(nrow(expr2),0.8*nrow(expr2),replace=F),]

tune.out2 <- tune.svm(if_health ~ .,data=training2,kernel='radial',decision.values=TRUE,
                     cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out2)
#bestmod<-tune.out$best.model
svmfit2 <- svm(if_health ~ .,data=training2,decision.values=TRUE,gamma=0.0001,cost=10)
summary(svmfit2)
test_pre2<-predict(svmfit2,testing2,decision.values=TRUE)
fitted2 <- attributes(test_pre2)$decision.values
roc_clu2 <- plot.roc(as.numeric(testing2$if_health),as.numeric(fitted2),
                    main=("ROC For GSE63898"), percent=TRUE, col="blue",print.auc=T)
table(true=testing2$if_health,pred=test_pre2)

##
gse63898_all_expr <- read.table("gse63898_all_expr.txt",sep="\t",header = T)
ran_gene2<-gse63898_all_expr[,sample((ncol(gse63898_all_expr)-1),334,replace = F)]
ran_data2<-cbind(ran_gene2,as.factor(gse63898_all_expr$if_health))
colnames(ran_data2)[ncol(ran_data2)]<-"if_health"

set.seed(12)
ran_training2 <- ran_data2[sample(nrow(ran_data2),0.8*nrow(ran_data2),replace=F),]
set.seed(12)
ran_testing2 <- ran_data2[-sample(nrow(ran_data2),0.8*nrow(ran_data2),replace=F),]

ran_tune.out2 <- tune.svm(if_health ~ .,data=ran_training2,kernel='radial',decision.values=TRUE,
                         cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
summary(ran_tune.out2)
#bestmod<-tune.out$best.model
ran_svmfit2 <- svm(if_health ~ .,data=ran_training2,decision.values=TRUE,gamma=0.0001,cost=10)
summary(ran_svmfit2)
ran_test_pre2<-predict(ran_svmfit2,ran_testing2,decision.values=TRUE)
ran_fitted2 <- attributes(ran_test_pre2)$decision.values
ran_roc_clu2 <- plot.roc(as.numeric(ran_testing2$if_health),as.numeric(ran_fitted2),
                        main=("ROC For GSE63898 random genes"), percent=TRUE, col="blue",print.auc=T)
table(true=ran_testing2$if_health,pred=ran_test_pre2)

#gse22058
##################
expr3<-data_22058
set.seed(125)
training3 <- expr3[sample(nrow(expr3),0.8*nrow(expr3),replace=F),]
set.seed(125)
testing3 <- expr3[-sample(nrow(expr3),0.8*nrow(expr3),replace=F),]

tune.out3 <- tune.svm(if_health ~ .,data=training3,kernel='radial',decision.values=TRUE,
                     cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
summary(tune.out3)
#bestmod<-tune.out$best.model
svmfit3 <- svm(if_health ~ .,data=training3,decision.values=TRUE,gamma=0.001,cost=0.1)
summary(svmfit3)
test_pre3<-predict(svmfit3,testing3,decision.values=TRUE)
fitted3 <- attributes(test_pre3)$decision.values
roc_clu3 <- plot.roc(as.numeric(testing3$if_health),as.numeric(fitted3),
                    main=("ROC For GSE22058"), percent=TRUE, col="blue",print.auc=T)
table(true=testing3$if_health,pred=test_pre3)

##
gse22058_all_expr <- read.table("gse22058_all_expr.txt",sep="\t",header = T)
ran_gene3<-gse22058_all_expr[,sample((ncol(gse22058_all_expr)-1),334,replace = F)]
ran_data3<-cbind(ran_gene3,as.factor(gse22058_all_expr$if_health))
colnames(ran_data3)[ncol(ran_data3)]<-"if_health"

set.seed(12)
ran_training3 <- ran_data3[sample(nrow(ran_data3),0.8*nrow(ran_data3),replace=F),]
set.seed(12)
ran_testing3 <- ran_data3[-sample(nrow(ran_data3),0.8*nrow(ran_data3),replace=F),]

ran_tune.out3 <- tune.svm(if_health ~ .,data=ran_training3,kernel='radial',decision.values=TRUE,
                         cost=c(0.1,0.5,1,5,10,20,100),gamma=c(0.0001,0.001,0.01,0.1))
summary(ran_tune.out3)
bestmod3<-ran_tune.out3$best.model
#ran_svmfit3 <- svm(if_health ~ .,data=ran_training3,decision.values=TRUE,gamma=0.0001,cost=5)
#summary(ran_svmfit3)
ran_test_pre3<-predict(bestmod3,ran_testing3,decision.values=TRUE)
ran_fitted3 <- attributes(ran_test_pre3)$decision.values
ran_roc_clu3 <- plot.roc(as.numeric(ran_testing3$if_health),as.numeric(ran_fitted3),
                        main=("ROC For GSE22058 random genes"), percent=TRUE, col="blue",print.auc=T)
table(true=ran_testing3$if_health,pred=ran_test_pre3)



#######
m<-ncol(expr)
svm<-NULL
for(i in 1:nrow(expr))
{
  label<-expr[-i,m]
  train_data<-expr[-i,1:m-1]
  pred_data<-expr[i,1:m-1]
  hcc_svm <- svm(train_data,label,type="C-classification",kernel="radial")
  pred <- predict(hcc_svm,pred_data,decision.values=TRUE)
  p<-attr(pred ,"decision.values")
  #sp<-sigmoid(p)
  svm<-rbind(svm,p)
}

svmf<-data.frame(svm)
pred2<-prediction(svmf,expr[,m])

b<-performance(pred2,measure="auc")
auc<-as.numeric(b@y.values)
auc

roc_clu<-plot.roc(as.numeric(expr[,m]),as.numeric(svm),
                  main=("ROC For test"),percent=TRUE, col="blue",print.auc=T)





#########
# fileName <- dir(pattern='GeneData_SelectbyRFE_ofCluster [0-9] .txt')
# auc9<-c()
# for(k in 1:length(fileName)){
  # data <- read.table(file = paste(fileName[k]),sep="\t",header = T)
  # # min_max<-function(x){ (x-min(x))/(max(x)-min(x))  }
  # # data_nor<-apply(data, 2, min_max)
  # 
  # mylabel<-c(rep(1,243),rep(0,268))
  # expr<-data.frame(data,if_health=as.factor(mylabel))
  # 
  # set.seed(12)
  # training <- expr[sample(nrow(expr),0.8*nrow(expr),replace=F),]
  # set.seed(12)
  # test <- expr[-sample(nrow(expr),0.8*nrow(expr),replace=F),]
  # # tune.out <- tune.svm(if_health ~ .,data=training,kernel='radial',
  # #                      cost=c(0.1,0.5,1,10,100),gamma=c(0.5,1,2,3,4))
  # # summary(tune.out)
  # # pred <- predict(tune.out$best.model,newx=dat[-train,])
  # # table(true=dat[-train,]$y,pred)
  # #,gamma=0.5,cost=1
  # svmfit <- svm(if_health ~ .,data=training,lernel='radial',cost=2,gamma=0.5,decision.values=TRUE)
  # fitted <- attributes(predict(svmfit,test,decision.values=TRUE))$decision.values
  # roc_clu <- plot.roc(as.numeric(test$if_health),as.numeric(fitted),
  #                     main=paste("ROC For cluster",k), percent=TRUE, col=col_roc[k],print.auc=T)
  # col_roc<-rainbow(9)
  # pdf(file=paste("ROC_Cluster",k,'.pdf',sep=''))
  # plot(roc_clu, main=paste("ROC For cluster",k), percent=TRUE, col=col_roc[k],print.auc=T)
  # dev.off()
  # auc1<-as.numeric( performance(prediction(fitted,test$if_health),measure="auc")@y.values)
  # auc9<-c(auc9,auc1)
# }

# hcc_svm <-svm(if_health ~ .,data = training,type="C-classification",kernel="radial")
# pre.forest<-predict(hcc_svm, test,decision.values=TRUE)
# table<-table(pre.forest,test$if_health)
# 
# p1<-attr(pre.forest ,"decision.values")
# sp1<-sigmoid(p1)
# perf1<-performance( prediction(sp1,test$if_health) ,"tpr","fpr" )
# auc1<-as.numeric( performance(prediction(sp1,test$if_health,measure="auc")@y.values) )
# auc1



