setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_feature\\tumor')
t_all<-read.table("tumor_all.txt",sep="\t",header = T)

tt<-t_all
row.names(tt)<-tt[,1]
t_all2<-tt[,-1]
#select the 
tumor_select<-t_all2[,c(1,3,4,5,6,7,11,12,13,14,15,16,19,20,23,24,28,29,31,33,36,40)]
write.table(tumor_select,"tumor_select.txt",quote=F,sep="\t")


