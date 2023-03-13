library(readxl)
library(dplyr)
library(data.table)
library(qpcR)

# DESeq2

# SNF2h-cohesin-NuRD complex
SNF2h_cohesin_NuRD=c("BAZ1A","CHD3","HDAC1","HDAC2","MBD2","MBD3","MTA1","MTA2","RAD21",
                     "RBBP4","RBBP7","SMARCA5","SMC1A","SMC3","STAG1","STAG2")

# ATRX_DAXX complex
ATRX_DAXX=c("ATRX", "DAXX")

# Cohesin_SA1/SA2 complex
Cohesin_SA1_SA2=c("RAD21","SMC1A","SMC3","STAG1","STAG2")

# SMAD4_SNO_SKI complex
SMAD4_SNO_SKI=c("SKI", "SKIL","SMAD4")

# IKBA_p50_p65 Nf(kappa)B complex
IKBA_p50_p65=c("NFKB1","NFKBIA","RELA")

variable_list=list(SNF2h_cohesin_NuRD,ATRX_DAXX,Cohesin_SA1_SA2,SMAD4_SNO_SKI, IKBA_p50_p65)
names(variable_list)=c("SNF2h_cohesin_NuRD","ATRX_DAXX","Cohesin_SA1_SA2","SMAD4_SNO_SKI","IKBA_p50_p65")


folders=c("gpath","plaq_n","nft","tangles","cogn_global_lv","cogng_random_slope")

for(v in 1:length(variable_list)){
  print(names(variable_list[v]))
  
  cogn=data.frame()
  
  for (j in 1:length(folders)){
    
    print(folders[j])
    
    variable=read.csv(paste0("G:/PFC_429_final_Spring2021/Ghada/ROSMAP_RNAseq/DESeq2/Results_z/RNAseq_",folders[j],".csv"))
    variable=variable[variable$Gene.name %in% variable_list[[v]],]
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
  
  #write.csv(cogn,file=paste0("G:/PFC_429_final_Spring2021/Ghada/CORUM/coreComplexes/bulk_RNA_DESeq2/",names(variable_list[v]),"_bulk_RNA_DESeq2.csv"))
  
}
