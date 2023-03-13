library(Seurat)
library(stringr)
library(dplyr)
setwd("D:/PFC_427_WS/Vanshika/muscat")
pbmc=readRDS("Vasc_integrated_batch_5000_predicted_final.rds")
all_genes=rownames(pbmc)
table=as.data.frame(all_genes)

setwd("D:/PFC_427_WS/Vanshika/muscat/Results/gpath")

filelist=list.files()

#upregulated genes

for (i in 1:length(filelist)){
  filename=filelist[[i]]
  name=str_remove_all(filelist[[i]], ".csv")
  print(filename)
  data=read.csv(filename)
  data_up=subset(data, logFC > 0)
  data_up=subset(data_up, p_adj.loc < 0.05)
  #data_down=subset(data, logFC < 0)
  #data_down=subset(data_down, p_adj.loc < 0.05)
  
  for (j in 1:nrow(table)){
    gene=table[j,"all_genes"]
    if (gene %in% data_up[,"gene"]){
      table[j,"additional_column"]<- 1
    } else {
      table[j,"additional_column"]<- 0
    }
  }
  names(table)[names(table) == "additional_column"] <- name
}
setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/gpath")
write.csv(table,file="consensus_table_upregulated.csv")

#downregulated genes
all_genes=rownames(pbmc)
table=as.data.frame(all_genes)

setwd("D:/PFC_427_WS/Vanshika/muscat/Results/gpath")

filelist=list.files()


for (i in 1:length(filelist)){
  filename=filelist[[i]]
  name=str_remove_all(filelist[[i]], ".csv")
  print(filename)
  data=read.csv(filename)
  data_down=subset(data, logFC < 0)
  data_down=subset(data_down, p_adj.loc < 0.05)
  
  for (j in 1:nrow(table)){
    gene=table[j,"all_genes"]
    if (gene %in% data_down[,"gene"]){
      table[j,"additional_column"]<- 1
    } else {
      table[j,"additional_column"]<- 0
    }
  }
  names(table)[names(table) == "additional_column"] <- name
}
setwd("D:/PFC_427_WS/Vanshika/muscat/DEG_consensus_signatures/gpath")
write.csv(table,file="consensus_table_downregulated.csv")