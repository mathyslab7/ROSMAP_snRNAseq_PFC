library(qpcR)

setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/cogn_global_lv/Exc")

load("consensus_gene_list_upreg.RData")

table=data.frame()
for(i in 1:length(consensus_gene_list)){
  df=as.data.frame(consensus_gene_list[i])
  if(length(table)>0) table=qpcR:::cbind.na(table,df)
  if(length(table)==0)table=df
}


write.csv(table, "consensus_gene_table_upregulated.csv")
