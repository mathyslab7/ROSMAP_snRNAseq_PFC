library(readxl)
library(dplyr)
library(data.table)
library(qpcR)

# gpath DGE - muscat

setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results")
folders=list.files()
folders=folders[!folders %like% ".RData"] 
folders=folders[folders %like% "corrected"]

for (j in 1:length(folders)){
  
  print(folders[j])
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/",folders[j]))
  folders_2=list.files()
  folders_2=folders_2[!folders_2 %like% ".RData"]
  
  for(k in 1:length(folders_2)){
    setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/",folders[j],"/",folders_2[k]))
    files=list.files()
    files=files[files %like% "Exc"]
    
    cogn=data.frame()
    for(f in 1:length(files)){
      celltype=read.csv(files[f])
      celltype=celltype[celltype$gene %in% c("HES4",
                                             "PDE10A",
                                             "RPH3A",
                                             "ST6GAL2",
                                             "UST"), ]
      celltype=celltype[order(celltype$gene),]
      
      celltype$logp=-log10(celltype$p_adj.loc)
      celltype$signed_logp=celltype$logp*(celltype$logFC/(abs(celltype$logFC)))
      celltype=celltype[,c("gene","signed_logp")]
      
      celltype_name=gsub(folders_2[k],"",files[f])
      celltype_name=gsub("_[^_]+$","",celltype_name)
      colnames(celltype)[2]=celltype_name
      
      if(length(cogn)>0) cogn=left_join(cogn,celltype,by="gene")
      if(length(cogn)==0) cogn=celltype
      
    }
    
    row.names(cogn)=cogn[,1]
    cogn=cogn[,-1]
    write.csv(cogn,file=paste0("G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/muscat/Association_scores/Exc/Exc_resilience_",folders[j],"_",folders_2[k],".csv"))
    
  }
}

