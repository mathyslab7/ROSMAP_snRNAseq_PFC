library(readxl)
library(dplyr)
library(data.table)
library(qpcR)

# gpath DGE - muscat

setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results")

#folders=list.files()
#folders=folders[!folders %like% ".RData"]
#folders=folders[grep("^cogn_|slope|CR|CDR",folders)]

folders=c("cogn_global_lv", "cogn_ep_lv", "cogn_po_lv", "cogn_ps_lv", "cogn_se_lv", "cogn_wo_lv")
folders=c("cogng_random_slope","cognep_random_slope", "cognpo_random_slope","cognps_random_slope", "cognse_random_slope", "cognwo_random_slope")
folders=c("gpath_CR_score","plaq_n_CR_score","nft_CR_score","tangles_CR_score")
folders=c("gpath_CDR_score","plaq_n_CDR_score","nft_CDR_score","tangles_CDR_score")


for (j in 1:length(folders)){
  
  print(folders[j])
  
  setwd(paste0("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results/",folders[j]))
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
    
    celltype_name=gsub(folders[j],"",files[f])
    celltype_name=gsub("_[^_]+$","",celltype_name)
    colnames(celltype)[2]=celltype_name
    
    if(length(cogn)>0) cogn=left_join(cogn,celltype,by="gene")
    if(length(cogn)==0) cogn=celltype
    
  }
  
  row.names(cogn)=cogn[,1]
  cogn=cogn[,-1]
  write.csv(cogn,file=paste0("F:/PFC_429_final_Spring2021/Ghada/Cognition/Cognitive Resilience/muscat/Association_scores/Exc/Exc_resilience_",folders[j],".csv"))
  
}

