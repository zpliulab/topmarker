setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_feature\\nontumor')

#spl
spl<-read.table("n_spl.txt",sep="\t",header = T)
tt<-spl
row.names(tt)<-NULL
splf<-data.frame(gene=row.names(spl),tt)

#python_10
nontumor_10<-read.table("nont_10.txt",sep="\t",header = F)
colnames(nontumor_10)<-c("gene","degree","degree_centrality","average_neighbor_degree","betweenness","eigenvector","cluster_coefficient","triangles","pagerank","load_centrality","harmonic_centrality")

my11<-merge(splf,nontumor_10,by="gene",all=FALSE)
write.table(my11,"n_my11.txt",quote = FALSE,sep = "\t")

#igraph
p_igraph<-read.table("n_igraph.txt",header=T,sep = "\t")

#sna
p_sna<-read.table("n_sna_3.txt",header=T,sep = "\t")

#centiserve
p_centiserve<-read.table("n_centiserve_15.txt",header=T,sep = "\t")

#braingraph
p_braingraph<-read.table("n_braingraph.txt",header=T,sep="\t")

#bind
p_r<-cbind(p_centiserve,p_igraph,p_braingraph,p_sna)
tt<-p_r
row.names(tt)<-NULL
p_rr<-cbind(gene=row.names(p_r),tt)

p_all<-merge(my11,p_rr,by="gene",all=FALSE)
write.table(p_all,"nontumor_all.txt",quote = FALSE,sep = "\t",row.names=F)
# 
# p_all[is.infinite(p_all[,25]),25]<-NA
# p_24_naomit<-na.omit(p_24)
# write.table(p_24_naomit,"tumor_24_naomit.txt",quote = FALSE,sep = "\t",row.names=F)







