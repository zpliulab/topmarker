# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 19:59:52 2018

@author: Wang Aki
"""

import os
os.chdir('C:\\Users\\Yanqiu Wang\\Desktop\\HCC\\3_network\\my_feature\\tumor') 

def mytwelve(infile,outfile):
    import networkx as nx
    open_file=open(infile,'r')
    myG=nx.Graph()     
    for line in open_file:
       head, tail = [str(x) for x in line.split()]  
       myG.add_edge(head,tail)
   
    degree = nx.degree(myG)
    d = nx.degree_centrality(myG)
    a = nx.average_neighbor_degree(myG)
    b = nx.betweenness_centrality(myG,normalized=False)
    e = nx.eigenvector_centrality(myG)
    lcc = nx.clustering(myG)
    t = nx.triangles(myG)
    pr = nx.pagerank(myG)
    #katz = nx.katz_centrality(myG,normalized=False)
    #infor = nx.information_centrality(myG)
    load = nx.load_centrality(myG,normalized=False)
    harmonic = nx.harmonic_centrality(myG)
    
    
    out_file=open(outfile,'w')
    for v in myG.nodes():
        out_file.write("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" %(v,degree[v],d[v],a[v],b[v],e[v],lcc[v],t[v],pr[v],load[v],harmonic[v]))
    
    out_file.close()
    open_file.close()    
   
mytwelve('t_cor_2gene.txt','t_10.txt')    
   
