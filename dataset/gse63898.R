setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_otherDATA2\\GSE63898')

#����gse63898���ݣ�ȡ�����б������ݺ�20����Ⱥ����880������ı������ݡ�
BigClu_gene_id<-read.table("BigClu_gene_id.txt",header=T,sep="\t")
gse_data<-read.table("GSE63898.gene.txt",sep="\t",header = T,fill=TRUE,quote="")
gse_BigClu<-merge(BigClu_gene_id,gse_data,by.x="GeneID",by.y="Entrez_ID",all=FALSE)
write.table(gse_BigClu,"gse63898_BigClu.txt",quote=F,sep="\t",row.names=F)#��880��BigClu�к��еĻ�����875���ཻ����

HCCDB13.sample<-read.table("HCCDB13.sample.txt",sep="\t",header = F,fill=TRUE,quote="")
sample<-t(HCCDB13.sample)
colnames(sample)<-sample[1,]
sam<-data.frame(sample[-1,])

t_sample<-sam[sam$TYPE=="HCC",]
n_sample<-sam[sam$TYPE=="Adjacent",]
t_sampleID<-gsub("-",".",t_sample$SAMPLE_ID)
n_sampleID<-gsub("-",".",n_sample$SAMPLE_ID)



######
#ȡ��875�����������
t_data<-gse_BigClu[,t_sampleID]
row.names(t_data)<-gse_BigClu$gene
t_gse_BigClu<-t(t_data)
n_data<-gse_BigClu[,n_sampleID]
row.names(n_data)<-gse_BigClu$gene
n_gse_BigClu<-t(n_data)

mylabel<-c(rep(1,dim(t_gse_BigClu)[1]),rep(0,dim(n_gse_BigClu)[1]))
gse63898_data<-rbind(n_gse_BigClu,t_gse_BigClu)
gse_expr<-data.frame(gse63898_data,if_health=as.factor(mylabel))
write.table(gse_expr,"gse63898_expr.txt",quote=F,sep="\t")


######
#ȡ�����л����������
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
gse63898_all_data<-rbind(n_gse_all,t_gse_all)
gse_all_expr<-data.frame(gse63898_all_data,if_health=as.factor(mylabel2))
write.table(gse_all_expr,"gse63898_all_expr.txt",quote=F,sep="\t")






