setwd('D:\\A_My_Data\\HCC\\1_data\\protein_validation')

homo<-read.table(file="All_Mammalia.gene_info",sep="\t",header=T,fill=TRUE,quote="")

colnames(homo)[1:3]<-c('tax_id','GeneID','Symbol')
homo2<-subset.data.frame(homo,select=c('tax_id','GeneID','Symbol'))
homo3<-homo2[homo2$tax_id == "10090",]
write.table(homo3,"mouse.txt",sep='\t',quote=F,row.names=F)

