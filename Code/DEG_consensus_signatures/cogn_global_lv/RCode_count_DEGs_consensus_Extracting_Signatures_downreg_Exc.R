setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/cogn_global_lv")


downregulated=read.csv("consensus_table_downregulated.csv")
rownames(downregulated)=downregulated[,2]
downregulated_Exc=downregulated[,c(9:22)]

downregulated_Exc$Count=rowSums(downregulated_Exc)
downregulated_Exc=cbind(Gene_Names=rownames(downregulated_Exc),downregulated_Exc)

downregulated_Exc2=downregulated_Exc[order(-downregulated_Exc$Count),]
downregulated_Exc2=subset(downregulated_Exc2,Count>1)

list_count=unique(downregulated_Exc2[,"Count"])

file_list=list.files(,pattern ="consensus_table")
name=gsub("consensus_table_","",file_list[1])
name=gsub(".csv","",name)

consensus_gene_list=list()

for(j in 1:length(list_count)){
  downregulated_Exc3=subset(downregulated_Exc2,Count>=list_count[[j]])
  gene_list=list(downregulated_Exc3[,"Gene_Names"])
  consensus_gene_list[j]=gene_list
  names(consensus_gene_list)[j]=paste0("consensus_",name,"_",list_count[j])
}
  
setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/cogn_global_lv/Exc")
save(consensus_gene_list,file="consensus_gene_list_downreg.RData")

#load("consensus_gene_list_downreg.RData")


