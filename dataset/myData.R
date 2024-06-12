setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\myData')

BigClu_gene_id<-read.table("BigClu_gene_id.txt",header=T,sep="\t")
gse_data<-read.table("GSE25097.gene.txt",sep="\t",header = T,fill=TRUE,quote="")
gse_BigClu<-merge(BigClu_gene_id,gse_data,by.x="GeneID",by.y="Entrez_ID",all=FALSE)
write.table(gse_BigClu,"gse25097_BigClu.txt",quote=F,sep="\t",row.names=F)#与880个BigClu中含有的基因有880个相交基因。

HCCDB3.sample<-read.table("HCCDB3.sample.txt",sep="\t",header = F,fill=TRUE,quote="")
sample<-t(HCCDB3.sample)
colnames(sample)<-sample[1,]
sam<-data.frame(sample[-1,])

t_sample<-sam[sam$TYPE=="HCC",]
n_sample<-sam[sam$TYPE=="Adjacent",]
t_sampleID<-gsub("-",".",t_sample$SAMPLE_ID)
n_sampleID<-gsub("-",".",n_sample$SAMPLE_ID)

t_data<-gse_BigClu[,t_sampleID]
row.names(t_data)<-gse_BigClu$gene
t_gse_BigClu<-t(t_data)
n_data<-gse_BigClu[,n_sampleID]
row.names(n_data)<-gse_BigClu$gene
n_gse_BigClu<-t(n_data)

mylabel<-c(rep(1,dim(t_gse_BigClu)[1]),rep(0,dim(n_gse_BigClu)[1]))
gse25097_data<-rbind(n_gse_BigClu,t_gse_BigClu)
gse_expr<-data.frame(gse25097_data,if_health=as.factor(mylabel))
write.table(gse_expr,"gse25097_expr.txt",quote=F,sep="\t")

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
gse25097_all_data<-rbind(n_gse_all,t_gse_all)
gse_all_expr<-data.frame(gse25097_all_data,if_health=as.factor(mylabel2))
write.table(gse_all_expr,"gse25097_all_expr.txt",quote=F,sep="\t")

