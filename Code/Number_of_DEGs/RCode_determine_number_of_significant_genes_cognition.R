library(Seurat)
library(muscat)
library(dplyr)

setwd("D:/PFC_427_WS/Vanshika/muscat/Results")

#Select variables
variable_list=list(
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
  "cogdx_AD_yes",
  "gpath_CR_score",
  "nft_CR_score",
  "plaq_n_CR_score",
  "tangles_CR_score",
  "gpath_CDR_score",
  "nft_CDR_score",
  "plaq_n_CDR_score",
  "tangles_CDR_score"
)

for (i in 1:length(variable_list)){
    if(grepl("cogdx_AD_yes",variable_list[[i]],fixed=TRUE)){
      variable="cogdx"
    } else {
      variable=variable_list[[i]]
    }  
    print(variable)
    path=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",variable_list[[i]])
    setwd(path)
    files=list.files(pattern="*.csv")
    for (j in 1:length(files)){
      pattern=paste0("_",variable,".csv")
      celltype=gsub(pattern,"",files[[j]])
      data=read.csv(files[[j]])
      data=subset(data,p_adj.loc < 0.05)
      nDEGs=nrow(data)
      nDEGs=as.data.frame(nDEGs)
      colnames(nDEGs)=celltype
      if (j==1){
        tempDF=nDEGs
      }
      if (j>1){
        tempDF=cbind(tempDF,nDEGs)
      }
    }
    rownames(tempDF)=variable_list[[i]]
    number_of_DEGs=t(tempDF)
    number_of_DEGs=cbind(cell_type_high_resolution=rownames(number_of_DEGs),number_of_DEGs)
    number_of_DEGs=as.data.frame(number_of_DEGs)
    if (i==1){
      results_table=number_of_DEGs
    }
    if (i>1){
      results_table=full_join(results_table,number_of_DEGs,by="cell_type_high_resolution")
    }
}

setwd("D:/PFC_427_WS/Vanshika/muscat/Number_of_significant_DEGs")
write.csv(results_table,file = "Number_of_significant_DEGs_cognition.csv")
