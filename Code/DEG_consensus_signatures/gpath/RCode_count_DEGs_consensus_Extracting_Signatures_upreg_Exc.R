setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/gpath")


upregulated=read.csv("consensus_table_upregulated.csv")
rownames(upregulated)=upregulated[,2]
upregulated_Exc=upregulated[,c(9:22)]

upregulated_Exc$Count=rowSums(upregulated_Exc)
upregulated_Exc=cbind(Gene_Names=rownames(upregulated_Exc),upregulated_Exc)

upregulated_Exc2=upregulated_Exc[order(-upregulated_Exc$Count),]
upregulated_Exc2=subset(upregulated_Exc2,Count>1)

list_count=unique(upregulated_Exc2[,"Count"])

file_list=list.files(,pattern ="consensus_table")
name=gsub("consensus_table_","",file_list[2])
name=gsub(".csv","",name)

consensus_gene_list=list()

for(j in 1:length(list_count)){
  upregulated_Exc3=subset(upregulated_Exc2,Count>=list_count[[j]])
  gene_list=list(upregulated_Exc3[,"Gene_Names"])
  consensus_gene_list[j]=gene_list
  names(consensus_gene_list)[j]=paste0("consensus_",name,"_",list_count[j])
}
  
setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/gpath/Exc")

save(consensus_gene_list,file="consensus_gene_list_upreg.RData")

#load("consensus_gene_list_upreg.RData")


