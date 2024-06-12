setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_feature\\nontumor')
n_all<-read.table("nontumor_all.txt",sep="\t",header = T)

tt<-n_all
row.names(tt)<-tt[,1]
n_all2<-tt[,-1]
#select the 
nontumor_select<-n_all2[,c(1,3,4,5,6,7,11,12,13,14,15,16,19,20,23,24,28,29,31,33,36,40)]
write.table(nontumor_select,"nontumor_select.txt",quote=F,sep="\t")


