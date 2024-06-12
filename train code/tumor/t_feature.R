setwd('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_feature\\tumor')
library(igraph)
library(centiserve)
library(brainGraph)
library(tidyr)
library(sna)
library(linkcomm)

ppi<-as.matrix(read.table("t_cor_2gene.txt",sep="\t",header = T))
pp<-graph.data.frame(ppi,directed=F)
ppa<-get.adjacency(pp,sparse=FALSE)

#braingraph
p_braingraph<-read.table("t_braingraph.txt",header=T,sep="\t")

#spl
spl<-read.table("t_spl.txt",sep="\t",header = T)
tt<-spl
row.names(tt)<-NULL
splf<-data.frame(gene=row.names(spl),tt)

#python_10
tumor_10<-read.table("t_10.txt",sep="\t",header = F)
colnames(tumor_10)<-c("gene","degree","degree_centrality","average_neighbor_degree","betweenness","eigenvector","cluster_coefficient","triangles","pagerank","load_centrality","harmonic_centrality")
my11<-merge(splf,tumor_10,by="gene",all=FALSE)
#write.table(my11,"my11.txt",quote = FALSE,sep = "\t",)


#igraph
i_ego_size<-data.frame(ego_size(pp))#
i_coreness<-data.frame(coreness(pp))
i_hub_score<-data.frame(hub_score(pp))
i_constraint<-data.frame(constraint(pp))
i_eccentricity<-data.frame(eccentricity(pp))
i_local_scan<-data.frame(local_scan(pp))
i_subgraph_centrality<-data.frame(subgraph_centrality(pp))
i_transitivity<-data.frame(transitivity(pp,type="local"))#

p_igraph<-cbind(i_ego_size,i_coreness,i_hub_score[,1],i_constraint,i_eccentricity,i_local_scan,i_subgraph_centrality,i_transitivity)
colnames(p_igraph)<-c("ego_size","coreness","hub_score","constraint","eccentricity","local_scan","subgraph_centrality","transitivity")
write.table(p_igraph,"t_igraph.txt",quote = FALSE,sep = "\t")


#sna
n_gilschmidt<-data.frame(gilschmidt(ppa,gmode="graph"))
n_prestige<-data.frame(prestige(ppa,gmode="graph"))
n_stresscent<-data.frame(stresscent(ppa,gmode="graph",cmode = "undirected"))

p_sna<-cbind(n_prestige,n_stresscent,n_gilschmidt)
colnames(p_sna)<-c("prestige","stresscent","gilschmidt")
write.table(p_sna,"t_sna_3.txt",quote = FALSE,sep = "\t")


#centiserve
c_markovcent<-data.frame(markovcent(pp))#
#c_averagedis<-data.frame(averagedis(pp))
#c_barycenter<-data.frame(barycenter(pp))
#c_closeness_cf<-data.frame(closeness.currentflow(pp))#
#c_closeness_freeman<-data.frame(closeness.freeman(pp))
#c_decay<-data.frame(decay(pp))
c_closeness_latora<-data.frame(closeness.latora(pp))
c_closeness_residual<-data.frame(closeness.residual(pp))
c_communitycent<-data.frame(communitycent(pp))
c_diffusion_degree<-data.frame(diffusion.degree(pp))
c_lobby<-data.frame(lobby(pp))
c_mnc<-data.frame(mnc(pp))
c_dmnc<-data.frame(dmnc(pp))
c_geokpath<-data.frame(geokpath(pp))
c_laplacian<-data.frame(laplacian(pp))
c_leverage<-data.frame(leverage(pp))
c_lincent<-data.frame(lincent(pp))
c_semilocal<-data.frame(semilocal(pp))
c_topocoefficient<-data.frame(topocoefficient(pp))
c_clusterrank<-data.frame(clusterrank(pp))#


p_centiserve<-cbind(c_markovcent,
                    c_closeness_latora,c_closeness_residual,c_communitycent,c_diffusion_degree,
                    c_lobby,c_mnc,c_dmnc,c_geokpath,c_laplacian,c_leverage,c_lincent,c_semilocal,
                    c_topocoefficient,c_clusterrank)
colnames(p_centiserve)<-c("markovcent","closeness_latora","closeness_residual",
                          "communitycent","diffusion_degree","lobby","mnc","dmnc",
                          "geokpath","laplacian","leverage","lincent","semilocal","topocoefficient","clusterrank")

write.table(p_centiserve,"t_centiserve_15.txt",quote = FALSE,sep = "\t",row.names = TRUE)


#bind
p_centiserve<-read.table("t_centiserve_8.txt",header=T,sep="\t")
p_igraph<-read.table("t_igraph.txt",header=T,sep="\t")


p_r<-cbind(p_centiserve,p_igraph,p_braingraph)
tt<-p_r
row.names(tt)<-NULL
p_rr<-cbind(gene=row.names(p_r),tt)

p_24<-merge(my9,p_rr,by="gene",all=FALSE)
p_23<-p_24[,-25]
write.table(p_23,"tumor_23.txt",quote = FALSE,sep = "\t",row.names=F)
p_24[is.infinite(p_24[,25]),25]<-NA
p_24_naomit<-na.omit(p_24)
write.table(p_24_naomit,"tumor_24_naomit.txt",quote = FALSE,sep = "\t",row.names=F)







