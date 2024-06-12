setwd('D:\\A_My_Data\\HCC\\3_network\\my_community\\my_comm_2')
library(igraph)
#community detection

ppi<-read.table("hcc_pw_gene_2.txt",sep="\t",header = T)
pp<-graph.data.frame(ppi,directed=F)

cfg<-cluster_fast_greedy(pp, merges = TRUE, modularity = TRUE, membership = TRUE)
# cl<-cluster_louvain(pp)
#找出两种分社群方法分出的社群之间的交叉基因数目
#有表格数据
#####
# overlap_gene_number<-matrix(0,nrow = 62,ncol = 65)
# cfg_number<-c()
# cl_number<-c()
# for (i in 1:62) {
#   nodes <- V(pp)[cfg$membership == i]
#   vertice1<-as.matrix(nodes)
#   v1<-data.frame(row.names(vertice1))
#   colnames(v1)<-"cfg"
#   cfg_number<-c(cfg_number,dim(v1)[1])
#   for (j in 1:65) {
#     nodes2 <- V(pp)[cl$membership == j]
#     vertice2<-as.matrix(nodes2)
#     v2<-data.frame(row.names(vertice2))
#     colnames(v2)<-"cl"
#     cl_number<-c(cl_number,dim(v2)[1])
#     overlap_gene<-merge(v1,v2,by.x = "cfg",by.y = "cl",all=FALSE)
#     overlap_gene_number[i,j]<-dim(overlap_gene)[1]
#   }
# }
# row.names(overlap_gene_number)<-cfg_number
# colnames(overlap_gene_number)<-cl_number[1:65]
# aa<-data.frame(cl_number[1:65])
# aa2<-cbind(aa,1:65)
# colnames(aa2)<-c("cl_number","col")
# ta<-aa2[order(aa2$cl_number,decreasing = TRUE),]
# ogn<-overlap_gene_number[,ta[,2]]
# write.table(overlap_gene_number,"overlap_gene_number.txt",quote = F,sep="\t")
# write.table(ogn,"ogn.txt",quote = F,sep="\t")



#####
#写出分出的社群各自包含的基因名文件。有58个社群，其中有20个基因数大于10个的社区。共971个基因。
member1<-as.matrix(membership(cfg))
gene1<-data.frame(cbind(row.names(member1),member1))
row.names(gene1)<-NULL
colnames(gene1)<-c("gene","cluster1")
# c62<-gene2[order(gene2$X2,decreasing = FALSE),]
comm_num<-data.frame()
for (i in 1:length(cfg)) {
  aa<-data.frame(gene1[gene1$cluster1==i,])
  comm_num<-rbind(comm_num,aa)
}
write.table(comm_num,"comm_num.txt",quote=F,sep="\t",row.names = F)

#33个hcc相关基因与社区基因有没有交叉。
homo_sapiens<-read.table(file="Homo_sapiens.gene_info",sep="\t",header=T,fill=TRUE,quote="")
homo<-subset.data.frame(homo_sapiens,select=c("GeneID","Symbol"))
comm_num_ID<-merge(comm_num,homo,by.x="gene",by.y="Symbol",all=FALSE)

disease33<-read.table("disease_33.txt",header=F,sep="\t")
dis_33gene<-data.frame(disease33$V2)
colnames(dis_33gene)<-"dis_geneid"
over_gene<-merge(dis_33gene,comm_num_ID,by.x="dis_geneid",by.y="GeneID",all=FALSE)
write.table(over_gene,"PWgene_33HCCgene_overlap.txt",quote=F,sep="\t",row.names = F)#2 overlap gene: CASP8,TERT


#写出每个社群中包含的基因
EachClu_num<-c()
for (k in 1:length(cfg)) {
  a<-data.frame(comm_num[comm_num$cluster1==k,1])
  colnames(a)<-"gene"
  tt<-dim(a)[1]
  EachClu_num<-c(EachClu_num,tt)
  write.table(a,file=paste("cluster",k,".txt"),sep="\t",quote=F,row.names=FALSE)
}
EachClu_numf<-data.frame(EachClu_num)
write.table(EachClu_numf,"EachClu_num.txt",quote=F,sep="\t",row.names = F)

#取出基因数大于10个的社群，放在一个数据框中。有20个社群。
tt<-1:length(cfg)
EachClu_numf2<-data.frame(cbind(EachClu_numf,tt))
mm<-EachClu_numf2[EachClu_numf2[,1]>=10,2]
bigclu_num<-EachClu_numf2[EachClu_numf2[,1]>=10,1]#20个基因数不小于10个的社群一共有880个基因。

BigClu<-data.frame()
for (j in 1:length(mm)) {
  aa<-data.frame(gene1[gene1$cluster1==mm[j],])
  BigClu<-rbind(BigClu,aa)
}
write.table(BigClu,"BigClu.txt",quote=F,sep="\t",row.names = F)


plot(cfg,pp,vertex.size =5,layout=layout.auto,vertex.label.cex=0.1,edge.color = grey(0.5))


#看筛选出的971基因与kegg_hcc_pathway中的基因有多少重合
library(tidyr)
# rfe_gene<-read.table("rfe_gene.txt",header = TRUE,sep="\t")
# rg<-as.data.frame(rfe_gene[,1])
# colnames(rg)<-"gene"

comm_node<-read.table("hcc_p_nodes_2.txt",header = TRUE,sep="\t")

hcc_pathway<-read.table("kegg_hcc_pathway.txt",header = TRUE,sep="\t")
colnames(hcc_pathway)<-c("ID","gene")
hp<-separate(hcc_pathway,col="gene",into=c("gene","v1"),sep = ";",remove=T,fill="right")
hp_gene<-data.frame(hp[,2])
colnames(hp_gene)<-"hp_gene"

share_gene<-merge(BigClu,hp_gene,by.x="gene",by.y="hp_gene",all=FALSE)
dim(share_gene)
share_gene
#167个kegg—hcc-gene，971个差异基因，只有13个交叉···
#share_gene: AKT2, E2F1, E2F2, FZD10, FZD5, LRP6, MTOR, PRKCG, TERT, TGFBR1, TXNRD1, WNT10B, WNT11

#20个社群与kegg—hcc-gene的交叉基因有7个: 
# gene    cluster1
# AKT2        2
# E2F1        14
# E2F2        2
# MTOR        2
# PRKCG       6
# TERT        2
# TGFBR1      4


######


bigclu<-read.table("BigClu.txt",header=T,sep="\t")
bigclu_gene<-bigclu[,1]
write.table(bigclu_gene,"bigclu_gene.txt",quote=F,sep="\t",row.names = F,col.names = F)







