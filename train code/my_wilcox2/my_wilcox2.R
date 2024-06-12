setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_net\\my_wilcox2')
library(igraph)
library(gtools)
#用网络特征来筛选差异基因，scale标准化，取交叉基因进行fold change，
#然后进行wilcoxon signed rank test，取p_value小于0.01的基因，
#取出的是疾病健康两个状态下网络特征差异比较大的基因。

t_p<-read.table("tumor_select.txt",header = T,sep = "\t")
t_p_scale<-data.frame(scale(t_p))

nont_p<-read.table("nontumor_select.txt",header = T,sep = "\t")
nont_p_scale<-data.frame(scale(nont_p))

#///////////////fc/////
t_gene<-data.frame(row.names(t_p_scale))
nont_gene<-data.frame(row.names(nont_p_scale))
colnames(t_gene)<-"gene"
colnames(nont_gene)<-"gene"
ol<-merge(t_gene,nont_gene,by.x="gene",by.y="gene",all=FALSE)

t_gene_p<-data.frame()
for (i in 1:nrow(ol)) {
  t<-row.names(t_p_scale)==ol[i,1] 
  a<-t_p_scale[t,]
  t_gene_p<-rbind(t_gene_p,a)
}
nont_gene_p<-data.frame()
for (i in 1:nrow(ol)) {
  t<-row.names(nont_p_scale)==ol[i,1] 
  a<-nont_p_scale[t,]
  nont_gene_p<-rbind(nont_gene_p,a)
}

t_mean<-data.frame(apply(t_gene_p,1,mean))
nt_mean<-data.frame(apply(nont_gene_p,1,mean))

# t<-c()
# for (i in 1:nrow(t_gene_p)) {
#   a=t_mean[i,1]/nt_mean[i,1]
#   if(a>=2 | a<=0.5){
#     t<-c(t,i)
#   }
# }
t<-c()
logFC<-c()
for (i in 1:nrow(t_gene_p)) {
  a = foldchange2logratio( foldchange(t_mean[i,1],nt_mean[i,1]), base=2) 
  logFC<-c(logFC,a)
  if( abs(a)>1 ){
    t<-c(t,i)
  }
}

logFC<-data.frame(logFC)

t_fc<-t_gene_p[t,]
nt_fc<-nont_gene_p[t,]


#////////////////WILCOXON RANK SUM////////////////////
# sigmoid<-function(x){  1/(1+exp(-x))  }
# t_sig<-data.frame(apply(t_fc,2,sigmoid))
# nont_sig<-data.frame(apply(nt_fc,2,sigmoid))

# tt_sig<-t(t_p_scale)
# tnont_sig<-t(nont_p_scale)

tt_fc<-t(t_fc)
tnt_fc<-t(nt_fc)


p_cor<-data.frame()
for(i in 1:nrow(t_fc)){
  a<-wilcox.test(tt_fc[,i],tnt_fc[,i],paired=TRUE)
  b<-data.frame(a$statistic,a$p.value)
  p_cor<-rbind(p_cor,b)
}

p_wilcox<-cbind(row.names(t_fc),p_cor)
colnames(p_wilcox)<-c("gene","w","p_value")
write.table(p_wilcox,"p_select_wilcox_2.txt",quote = FALSE,sep = "\t")

pw_order<-p_wilcox[order(p_wilcox$p_value,decreasing = FALSE),]
pw_min<-pw_order[pw_order$p_value<=0.01,]
write.table(pw_order,"pw_select_order_2.txt",quote = FALSE,sep = "\t")
write.table(pw_min,"pw_select_min_2.txt",quote = FALSE,sep = "\t")

##
tt<-rep(1,nrow(pw_min))
pw_Gene<-cbind(data.frame(pw_min$gene),tt)
colnames(pw_Gene)<-c('gene',"cluster")
write.table(pw_Gene,"pw_select_min_Gene_2.txt",quote = FALSE,sep = "\t",row.names = F)


#//////////////////////////////////////////////////////
gene_pw<-data.frame(pw_min[,1])
colnames(gene_pw)<-"gene"

t_2gene<-read.table("t_cor_2gene.txt",header = F,sep="\t")
nont_2gene<-read.table("nont_cor_2gene.txt",header = F,sep="\t")
colnames(t_2gene)<-c("gene1","gene2")
colnames(nont_2gene)<-c("gene1","gene2")

t1<-merge(gene_pw,t_2gene,by.x="gene",by.y="gene1",all=FALSE)
t2<-merge(gene_pw,t1,by.x="gene",by.y="gene2",all=FALSE)
tn1<-as.data.frame(t2[,1])
tn2<-as.data.frame(t2[,2])
colnames(tn1)<-"gene"
colnames(tn2)<-"gene"
tn<-rbind(tn1,tn2)
tn_u<-unique(tn)
print(sprintf("gene_pw_min: %d",nrow(gene_pw)))
print(sprintf("tumor nodes: %d",nrow(tn_u)))
print(sprintf("tumor edges: %d",nrow(t2)))

nt1<-merge(gene_pw,nont_2gene,by.x="gene",by.y="gene1",all=FALSE)
nt2<-merge(gene_pw,nt1,by.x="gene",by.y="gene2",all=FALSE)
ntn1<-as.data.frame(nt2[,1])
ntn2<-as.data.frame(nt2[,2])
colnames(ntn1)<-"gene"
colnames(ntn2)<-"gene"
ntn<-rbind(ntn1,ntn2)
ntn_u<-unique(ntn)
print(sprintf("nontumor nodes: %d",nrow(ntn_u)))
print(sprintf("nontumor edges: %d",nrow(nt2)))

u1<-as.data.frame(tn_u)
u2<-as.data.frame(ntn_u)
colnames(u1)<-"gene"
colnames(u2)<-"gene"
u_t_n<-merge(u1,u2,all=FALSE)#tumor nontumor overlap
print(sprintf("tumor and nontumor overlap ondes: %d",nrow(u_t_n)))
write.table(u_t_n,"u_t_n_2.txt",quote = FALSE,sep = "\t",row.names = F)

write.table(t2,"t_pw_gene_2.txt",quote = FALSE,sep = "\t",row.names = F)
write.table(nt2,"nont_pw_gene_2.txt",quote = FALSE,sep = "\t",row.names = F)
write.table(tn_u,"t_pw_nodes_2.txt",quote = FALSE,sep = "\t",row.names = F)
write.table(ntn_u,"nont_pw_nodes_2.txt",quote = FALSE,sep = "\t",row.names = F)


hcc<-read.table("hcc_genepairs.txt",header = T,sep="\t")
h1<-merge(gene_pw,hcc,by.x="gene",by.y="Gene1",all=FALSE)
h2<-merge(gene_pw,h1,by.x="gene",by.y="Gene2",all=FALSE)
write.table(h2,"hcc_pw_gene_2.txt",quote = FALSE,sep = "\t",row.names = F)
hn1<-as.data.frame(h2[,1])
hn2<-as.data.frame(h2[,2])
colnames(hn1)<-"gene"
colnames(hn2)<-"gene"
hn<-rbind(hn1,hn2)
hn_u<-unique(hn)
write.table(hn_u,"hcc_p_nodes_2.txt",quote = FALSE,sep = "\t",row.names = F)
print(sprintf("hcc nodes: %d",nrow(hn_u)))
print(sprintf("hcc edges: %d",nrow(h2)))






