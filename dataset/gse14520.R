setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\GSE14520')

#整理gse14520数据，取出所有表达数据和20个社群共有880个基因的表达数据。
BigClu_gene_id<-read.table("BigClu_gene_id.txt",header=T,sep="\t")
gse_data<-read.table("GSE14520-GPL3921.gene.txt",sep="\t",header = T,fill=TRUE,quote="")
gse_BigClu<-merge(BigClu_gene_id,gse_data,by.x="GeneID",by.y="Entrez_ID",all=FALSE)
write.table(gse_BigClu,"gse14520_BigClu.txt",quote=F,sep="\t",row.names=F)#与880个BigClu中含有的基因有770个相交基因。

HCCDB6.sample<-read.table("HCCDB6.sample.txt",sep="\t",header = F,fill=TRUE,quote="")
sample<-t(HCCDB6.sample)
colnames(sample)<-sample[1,]
sam<-data.frame(sample[-1,])

t_sample<-sam[sam$TYPE=="HCC",]
n_sample<-sam[sam$TYPE=="Adjacent",]
t_sampleID<-gsub("-",".",t_sample$SAMPLE_ID)
n_sampleID<-gsub("-",".",n_sample$SAMPLE_ID)
######
#取出770基因表达数据
t_data<-gse_BigClu[,t_sampleID]
row.names(t_data)<-gse_BigClu$gene
t_gse_BigClu<-t(t_data)
n_data<-gse_BigClu[,n_sampleID]
row.names(n_data)<-gse_BigClu$gene
n_gse_BigClu<-t(n_data)
mylabel<-c(rep(1,dim(t_gse_BigClu)[1]),rep(0,dim(n_gse_BigClu)[1]))
gse14520_data<-rbind(n_gse_BigClu,t_gse_BigClu)
gse_expr<-data.frame(gse14520_data,if_health=as.factor(mylabel))
write.table(gse_expr,"gse14520_expr.txt",quote=F,sep="\t")

######
#取出所有基因表达数据
homo<-read.table(file="Homo_sapiens.gene_info",sep="\t",header=T,fill=TRUE,quote="")
homo2<-subset.data.frame(homo,select=c("GeneID","Symbol"))
gse_homo<-merge(homo2,gse_data,by.x="GeneID",by.y="Entrez_ID",all=FALSE)

t_alldata<-gse_homo[,t_sampleID]
row.names(t_alldata)<-gse_homo$Symbol.x
t_gse_all<-t(t_alldata)
n_alldata<-gse_homo[,n_sampleID]
row.names(n_alldata)<-gse_homo$Symbol.x
n_gse_all<-t(n_alldata)

mylabel2<-c(rep(1,dim(t_gse_all)[1]),rep(0,dim(n_gse_all)[1]))
gse14520_all_data<-rbind(n_gse_all,t_gse_all)
gse_all_expr<-data.frame(gse14520_all_data,if_health=as.factor(mylabel2))
write.table(gse_all_expr,"gse14520_all_expr.txt",quote=F,sep="\t")







