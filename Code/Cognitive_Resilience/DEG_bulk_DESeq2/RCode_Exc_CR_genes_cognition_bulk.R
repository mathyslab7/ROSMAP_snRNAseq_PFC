library(readxl)
library(dplyr)
library(data.table)
library(qpcR)

# DESeq2

setwd("G:/PFC_429_final_Spring2021/Ghada/ROSMAP_RNAseq/DESeq2/Results_z")

folders=c("cogn_global_lv", "cogn_ep_lv", "cogn_po_lv", "cogn_ps_lv", "cogn_se_lv", "cogn_wo_lv")
folders=c("cogng_random_slope","cognep_random_slope", "cognpo_random_slope","cognps_random_slope", "cognse_random_slope", "cognwo_random_slope")
folders=c("gpath_CR_score","plaq_n_CR_score","nft_CR_score","tangles_CR_score")
folders=c("gpath_CDR_score","plaq_n_CDR_score","nft_CDR_score","tangles_CDR_score")


cogn=data.frame()

for (j in 1:length(folders)){
  
  print(folders[j])
  
  variable=read.csv(paste0("G:/PFC_429_final_Spring2021/Ghada/ROSMAP_RNAseq/DESeq2/Results_z/RNAseq_",folders[j],".csv"))
  variable=variable[variable$Gene.name %in% c("HES4",
                                              "PDE10A",
                                              "RPH3A",
                                              "ST6GAL2",
                                              "UST"
  ), ]
  variable=variable[order(variable$Gene.name),]
  
  variable$logp=-log10(variable$padj)
  variable$signed_logp=variable$logp*(variable$log2FoldChange/(abs(variable$log2FoldChange)))
  variable=variable[,c("Gene.name","signed_logp")]
  
  variable_name=folders[j]
  colnames(variable)[2]=variable_name
  
  if(length(cogn)>0) cogn=full_join(cogn,variable,by="Gene.name")
  if(length(cogn)==0) cogn=variable
  
}

row.names(cogn)=cogn[,1]
cogn=cogn[,-1]

write.csv(cogn,file="G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/edgeR_DESeq2/Association_scores//Exc/bulk_resilience_CDR_DESeq2.csv")

