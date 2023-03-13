library(stringr)
library(qpcR)
library(dplyr)
library(data.table)


celltypes=read.csv("G:/PFC_429_final_Spring2021/Ghada/Celltypes Annotations and Numbers/celltypes_order_muscat.csv")
celltypes=unlist(celltypes$celltype)
celltypes=gsub("/","_",celltypes)
celltypes=gsub("\\s*\\([^\\)]+\\)","",as.character(celltypes))
#celltypes=celltypes[1]
celltypes=celltypes[celltypes %like% "Exc"]

table_all=data.frame()

for(cell in 1:length(celltypes)){
  
  print(celltypes[cell])
  
  setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results")
  folders=list.files()
  folders=folders[!folders %like% ".RData"] 
  folders=folders[folders %like% "CR|CDR"]
  
  var_list=list()
  for (i in 1:length(folders)){
    
    setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/", folders[i]))
    files=list.files()
    
    files=files[files %like% celltypes[cell]]
    filename=gsub(".csv","",files)
    data=read.csv(files)
    
    data_up=subset(data, logFC > 0)
    data_up=subset(data_up, p_adj.loc < 0.05)
    
    data_up=data_up[,"gene"]
    name=folders[i]
    
    var_list[[name]]=data_up
    
  }
  
  CR_genes_up=unique(unlist(var_list))

#####################

  
  var_list_all=var_list
  genes_all=unique(CR_genes_up)
  
  
  if(length(genes_all>0)) {
  
  # counts 
  
  table=data.frame()
  for(g in 1:length(genes_all)){
    name=genes_all[g]
    count=0
    
    table_2=data.frame()
    for(v in 1:length(var_list_all)){
      if(genes_all[g] %in% var_list_all[[v]]){
        
        if(count>0) count=count+1
        if(count==0) count=1
      }
    }
    table_2=cbind(name,count)
    
    if(nrow(table)>0) table=rbind(table,table_2)
    if(nrow(table)==0) table=table_2 
  }
  
  
  table=as.data.frame(table)
  table$count=as.numeric(table$count)
  
  name=celltypes[cell]
  name=gsub(" ","_",name)
  name=gsub("-","_",name)
  name=gsub("/","_",name)
  
  setwd("G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/muscat/Summary/Gene_Scoring/Upreg")
  write.csv(table, file=paste0(name,"_DEG_upreg_resilience_genes_counts_all_conditions.csv"))
  
  table=table[table$count==8,]
  table=table[,1,drop=FALSE]
  colnames(table)=celltypes[cell]
  
  if(nrow(table)>=1){
  if(length(table_all)>0) table_all=qpcR:::cbind.na(table_all,table)
  if(length(table_all)==0) table_all=table
  }
  
  }
}


setwd("G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/muscat/Summary/Gene_Scoring/Upreg")
write.csv(table_all, file=paste0("Summary_DEG_upreg_resilience_genes_counts_CR_CDR.csv"))
  
  
  



