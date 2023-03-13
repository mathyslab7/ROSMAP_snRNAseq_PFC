setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results")
library(Seurat)
library(muscat)
library(dplyr)

#Select variables
variable_list=list(
  "amyloid",
  "gpath",
  "gpath_3neocort",
  "nft",
  "nft_mf",
  "tangles",
  "plaq_n",
  "plaq_n_mf",
  "plaq_d",
  "plaq_d_mf",
  "Apoe_e4",
  "arteriol_scler",
  "bradysc_lv",
  "caa_4gp",
  "cancer_bl",
  "ci_num2_gct",
  "ci_num2_gtt",
  "ci_num2_mct",
  "ci_num2_mtt",
  "cogdx_stroke",
  "cvda_4gp2",
  "diabetes_sr_rx_bl",
  "dlbdx",
  "dxpark",
  "gaitsc_lv",
  "headinjrloc_bl",
  "heart_bl",
  "hypertension_bl",
  "msex",
  "parksc_lv",
  "stroke_bl",
  "tdp_st4",
  "cogdx",
  "cognition",
  "cogn_global_lv",
  "cogn_ep_lv",
  "cogn_po_lv",
  "cogn_ps_lv",
  "cogn_se_lv",
  "cogn_wo_lv",
  "cogng_random_slope",
  "cognep_random_slope",
  "cognpo_random_slope",
  "cognps_random_slope",
  "cognse_random_slope",
  "cognwo_random_slope",
  "gpath_CR_score",
  "nft_CR_score",
  "plaq_n_CR_score",
  "cogdx_AD_yes",
 "tangles_CR_score",
 "chd_cogact_freq",
 "lifetime_cogact_freq_bl",
 "ma_adult_cogact_freq",
 "phys5itemsum_bl",
 "phys5itemsum_lv",
 "soc_net_bl",
 "social_isolation_avg",
 "social_isolation_lv",
 "ya_adult_cogact_freq"
)

variable=variable_list[[1]]

pattern=paste0("*",variable,".RData")
list.filenames=list.files(pattern=pattern)

#number_of_DEGs=data.frame()

load(list.filenames[1])
# access results table for 1st comparison
tbl <- res$table[[1]]
for (i in 1:length(tbl)) {
  tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
  tbl[[i]]=subset(tbl[[i]],logFC < 0)
  tbl[[i]]=nrow(tbl[[i]])
}
number_of_DEGs<-as.data.frame(tbl)

for (i in 2:length(list.filenames)){
  load(list.filenames[i])
  # access results table for 1st comparison
  tbl <- res$table[[1]]
  for (i in 1:length(tbl)) {
    tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
    tbl[[i]]=subset(tbl[[i]],logFC < 0)
    tbl[[i]]=nrow(tbl[[i]])
  }
  DF<-as.data.frame(tbl)
  number_of_DEGs=cbind(number_of_DEGs,DF)
}
rownames(number_of_DEGs)<-variable_list[[1]]
number_of_DEGs=t(number_of_DEGs)
number_of_DEGs=cbind(cell_type_high_resolution=rownames(number_of_DEGs),number_of_DEGs)
number_of_DEGs=as.data.frame(number_of_DEGs)


for (j in 2:length(variable_list)) {
  
  #Select variable (here cogdx Workspaces)
  variable=variable_list[[j]]
  
  pattern=paste0("*",variable,".RData")
  list.filenames=list.files(pattern=pattern)
  
  #number_of_DEGs=data.frame()
  
  load(list.filenames[1])
  # access results table for 1st comparison
  tbl <- res$table[[1]]
  for (i in 1:length(tbl)) {
    tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
    tbl[[i]]=subset(tbl[[i]],logFC < 0)
    tbl[[i]]=nrow(tbl[[i]])
  }
  number_of_DEGs2<-as.data.frame(tbl)
  
  for (i in 2:length(list.filenames)){
    load(list.filenames[i])
    # access results table for 1st comparison
    tbl <- res$table[[1]]
    for (i in 1:length(tbl)) {
      tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
      tbl[[i]]=subset(tbl[[i]],logFC < 0)
      tbl[[i]]=nrow(tbl[[i]])
    }
    DF<-as.data.frame(tbl)
    number_of_DEGs2=cbind(number_of_DEGs2,DF)
  }
  rownames(number_of_DEGs2)<-variable_list[[j]]
  number_of_DEGs2=t(number_of_DEGs2)
  number_of_DEGs2=cbind(cell_type_high_resolution=rownames(number_of_DEGs2),number_of_DEGs2)
  number_of_DEGs2=as.data.frame(number_of_DEGs2)
  
  number_of_DEGs=full_join(number_of_DEGs,number_of_DEGs2,by="cell_type_high_resolution")
}


setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Number_of_significant_DEGs")
write.csv(number_of_DEGs,file = "DGE_Comprehensive_number_of_significant_genes_down.csv")
