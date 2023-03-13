library(stringr)
library(qpcR)
setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")
filelist=list.files()

table=data.frame()
for (i in 1:length(filelist)){
  table=data.frame()
  filename=filelist[[i]]
  print(filename)
  setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")
  data=read.csv(filename)
  data_up=subset(data, logFC > 0)
  data_up=subset(data_up, p_adj.loc < 0.05)
  data_down=subset(data, logFC < 0)
  data_down=subset(data_down, p_adj.loc < 0.05)
  
  name=str_remove_all(filelist[[i]], ".csv")
  name_up=paste0(name,"_up")
  name_down=paste0(name,"_down")
  
  up_genes=data_up[,"gene",drop=FALSE]
  down_genes=data_down[,"gene",drop=FALSE]
  names(up_genes)[names(up_genes) == "gene"] <- name_up
  names(down_genes)[names(down_genes) == "gene"] <- name_down
  
  table=up_genes
  table <- qpcR:::cbind.na(table, down_genes)
  
  setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/GO_Enrichment/gpath")
  output_name=paste0(name,".csv")
  write.csv(table,file=output_name,row.names=FALSE)
}