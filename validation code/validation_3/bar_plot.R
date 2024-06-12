setwd('D:\\A_My_Data\\HCC\\1_data\\validation_3')
library(reshape2)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

#mydata <- melt(data1,id.vars="Conpany",variable.name="Year",value.name="Sale")

#top_clu<-read.table("t_test_AUC.txt",header=T,sep='\t')
# random33<-read.table("mean_auc_random33.txt",header=T,sep='\t')
# mean_auc_p_1_1000<-read.table("mean_auc_p_1_1000.txt",header=T,sep='\t')
# mean_auc_p_1001_2000<-read.table("mean_auc_p_1001_2000.txt",header=T,sep='\t')
# mean_auc_p_2001_3000<-read.table("mean_auc_p_2001_3000.txt",header=T,sep='\t')
# mean_auc_p_3001_4000<-read.table("mean_auc_p_3001_4000.txt",header=T,sep='\t')

top_clu<-read.table("limma_top33_AUC.txt",header=T,sep='\t')
random33<-read.table("limma_mean_auc_random33.txt",header=T,sep='\t')
mean_auc_p_1_1000<-read.table("limma_mean_auc_p_1_1000.txt",header=T,sep='\t')
mean_auc_p_1001_2000<-read.table("limma_mean_auc_p_1001_2000.txt",header=T,sep='\t')
mean_auc_p_2001_3000<-read.table("limma_mean_auc_p_2001_3000.txt",header=T,sep='\t')
mean_auc_p_3001_4000<-read.table("limma_mean_auc_p_3001_4000.txt",header=T,sep='\t')


all_auc<-cbind(random33,mean_auc_p_1_1000,mean_auc_p_1001_2000,mean_auc_p_2001_3000,mean_auc_p_3001_4000)
colnames(all_auc)<-c('random33','mean_auc_p_1_1000','mean_auc_p_1001_2000','mean_auc_p_2001_3000','mean_auc_p_3001_4000')
row.names(all_auc)<-c('GSE14520_auc','GSE22058_auc','GSE63898_auc','GSE64041_auc','GSE45436_auc','TCGA LIHC_auc')

mod3<-data.frame(top_clu$clu3_AUC)
colnames(mod3)<-"x"
top33<-data.frame(top_clu$t_test_AUC)
colnames(top33)<-"x"


t0<-cbind(mod3,scale="DNP 33 genes",data=row.names(all_auc))
t1<-cbind(top33,scale="DEGs top33 genes",data=row.names(all_auc))
t2<-cbind(mean_auc_p_1_1000,scale="P-value rank 1-1000",data=row.names(all_auc))
t3<-cbind(mean_auc_p_1001_2000,scale="P-value rank 1001-2000",data=row.names(all_auc))
t4<-cbind(mean_auc_p_2001_3000,scale="P-value rank 2001-3000",data=row.names(all_auc))
t5<-cbind(mean_auc_p_3001_4000,scale="P-value rank 3001-4000",data=row.names(all_auc))
t6<-cbind(random33,scale="Random33 genes",data=row.names(all_auc))


mydata<-rbind(t6,t5,t4,t3,t2,t1,t0)
colnames(mydata)[1]<-"AUC"

mypalette2<-brewer.pal(6,"Set3")
mypalette2<-c( "#FFFFB3", "#FDB462", "#FB8072" ,"#BEBADA","#80B1D3","#8DD3C7" )
mypalette2<-c("#8DD3C7" , "lightyellow","#FDB462", "#80B1D3", "#FB8072", "#BEBADA" )
mypalette2<-c( "lightcoral","lightsalmon" , "lemonchiffon","darkseagreen","lightseagreen", "lightblue" )


ggplot(mydata,aes(scale,AUC,fill=data))+
  geom_bar(stat="identity",position="dodge")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = mypalette2)+
  coord_flip()
  


