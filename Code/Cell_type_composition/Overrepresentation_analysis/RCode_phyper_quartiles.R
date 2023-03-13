library(dplyr)
setwd("D:/PFC_427_WS/Vanshika")
celltypes=read.csv("celltype_order_Mic_Ast.csv")
celltypes2=celltypes[!(celltypes$celltype=="Ast" | celltypes$celltype=="Mic" | celltypes$celltype=="Exc RELN CHD7" | celltypes$celltype=="Exc NRGN"),]

setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects")
data=read.csv("PFC427_meta.csv")

#Exclude Exc RELN CHD7 and Exc NRGN
data2<-data[!(data$cell_type_high_resolution=="Exc RELN CHD7" | data$cell_type_high_resolution=="Exc NRGN"),]

setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Confirmation_Cell_numbers")
metadata=read.csv("Combined_metadata_v2_quartile_median.csv")

select=c("projid",
  "cogn_global_lv_quartile",      
  "cogn_ep_lv_quartile",
  "cogn_po_lv_quartile",
  "cogn_ps_lv_quartile",                      
  "cogn_se_lv_quartile",
  "cogn_wo_lv_quartile",
  "cogng_random_slope_quartile",
  "cognep_random_slope_quartile",
  "cognpo_random_slope_quartile",
  "cognps_random_slope_quartile",
  "cognse_random_slope_quartile",
  "cognwo_random_slope_quartile",
  "gpath_quartile",                           
  "gpath_3neocort_quartile",
  "amyloid_quartile",
  "plaq_d_quartile",
  "plaq_d_mf_quartile",
  "plaq_n_quartile",                          
  "plaq_n_mf_quartile",
  "nft_quartile",
  "nft_mf_quartile",
  "tangles_quartile",
  "gpath_CDR_score_quartile",
  "plaq_n_CDR_score_quartile",               
  "nft_CDR_score_quartile",
  "tangles_CDR_score_quartile",
  "gpath_CR_score_quartile",
  "plaq_n_CR_score_quartile",
  "nft_CR_score_quartile",
  "tangles_CR_score_quartile",
  
  "cogn_global_lv_median",      
  "cogn_ep_lv_median",
  "cogn_po_lv_median",
  "cogn_ps_lv_median",                      
  "cogn_se_lv_median",
  "cogn_wo_lv_median",
  "cogng_random_slope_median",
  "cognep_random_slope_median",
  "cognpo_random_slope_median",
  "cognps_random_slope_median",
  "cognse_random_slope_median",
  "cognwo_random_slope_median",
  "gpath_median",                           
  "gpath_3neocort_median",
  "amyloid_median",
  "plaq_d_median",
  "plaq_d_mf_median",
  "plaq_n_median",                          
  "plaq_n_mf_median",
  "nft_median",
  "nft_mf_median",
  "tangles_median",
  "gpath_CDR_score_median",
  "plaq_n_CDR_score_median",               
  "nft_CDR_score_median",
  "tangles_CDR_score_median",
  "gpath_CR_score_median",
  "plaq_n_CR_score_median",
  "nft_CR_score_median",
  "tangles_CR_score_median"
  )

metadata2=metadata[select]

data_table=left_join(data2,metadata2,by="projid")

variable_list=c("cogn_global_lv_quartile",
                "cogn_ep_lv_quartile",
                       "cogn_po_lv_quartile",
                       "cogn_ps_lv_quartile",                      
                       "cogn_se_lv_quartile",
                       "cogn_wo_lv_quartile",
                        "cogng_random_slope_quartile",
                       "cognep_random_slope_quartile",
                       "cognpo_random_slope_quartile",
                       "cognps_random_slope_quartile",
                       "cognse_random_slope_quartile",
                       "cognwo_random_slope_quartile",
                       "gpath_quartile",                           
                       "gpath_3neocort_quartile",
                       "amyloid_quartile",
                       "plaq_d_quartile",
                       "plaq_d_mf_quartile",
                       "plaq_n_quartile",                          
                       "plaq_n_mf_quartile",
                       "nft_quartile",
                       "nft_mf_quartile",
                       "tangles_quartile",
                       "gpath_CDR_score_quartile",
                       "plaq_n_CDR_score_quartile",               
                       "nft_CDR_score_quartile",
                       "tangles_CDR_score_quartile",
                       "gpath_CR_score_quartile",
                       "plaq_n_CR_score_quartile",
                       "nft_CR_score_quartile",
                       "tangles_CR_score_quartile")

for (i in 1:length(variable_list)){
  variable=variable_list[i]
  selected_columns=c("cell_type_high_resolution",variable)
  data_table_temp=data_table[selected_columns]
  data_table_temp=na.omit(data_table_temp)
  total_cell_number=nrow(data_table_temp)
  
  total_first=sum(na.omit(data_table_temp[variable] == "first_quarter"))
  total_second=sum(na.omit(data_table_temp[variable] == "second_quarter"))
  total_third=sum(na.omit(data_table_temp[variable] == "third_quarter"))
  total_fourth=sum(na.omit(data_table_temp[variable] == "fourth_quarter"))
  
  for (j in 1:length(celltypes2)){
    cell_type=celltypes2[j]
    ncells_celltype=sum(na.omit(data_table_temp["cell_type_high_resolution"] == cell_type))
    celltype_df=subset(data_table_temp,cell_type_high_resolution==cell_type)
    ct_first=sum(na.omit(celltype_df[variable] == "first_quarter"))
    ct_second=sum(na.omit(celltype_df[variable] == "second_quarter"))
    ct_third=sum(na.omit(celltype_df[variable] == "third_quarter"))
    ct_fourth=sum(na.omit(celltype_df[variable] == "fourth_quarter"))
    
    first_result=phyper((ct_first-1),total_first,(total_cell_number-total_first),ncells_celltype,lower.tail= FALSE)
    second_result=phyper((ct_second-1),total_second,(total_cell_number-total_second),ncells_celltype,lower.tail= FALSE)
    third_result=phyper((ct_third-1),total_third,(total_cell_number-total_third),ncells_celltype,lower.tail= FALSE)
    fourth_result=phyper((ct_fourth-1),total_fourth,(total_cell_number-total_fourth),ncells_celltype,lower.tail= FALSE)
    
    first_result=as.data.frame(first_result)
    second_result=as.data.frame(second_result)
    third_result=as.data.frame(third_result)
    fourth_result=as.data.frame(fourth_result)
    
    rownames(first_result)=cell_type
    rownames(second_result)=cell_type
    rownames(third_result)=cell_type
    rownames(fourth_result)=cell_type
    
    colname_first=paste0(variable,"_first_quarter")
    colname_second=paste0(variable,"_second_quarter")
    colname_third=paste0(variable,"_third_quarter")
    colname_fourth=paste0(variable,"_fourth_quarter")
    
    colnames(first_result)=colname_first
    colnames(second_result)=colname_second
    colnames(third_result)=colname_third
    colnames(fourth_result)=colname_fourth
    
    temp_combined=cbind(first_result,second_result,third_result,fourth_result)
    temp_combined=cbind(celltype=rownames( temp_combined), temp_combined)
    if(j==1){
      temp_results=temp_combined
    }
    if(j>1){
      temp_results=rbind(temp_results,temp_combined)
    }
  }
  if (i==1){
    results_table=temp_results
  }
  if (i>1){
    results_table=full_join(results_table,temp_results,by="celltype")
  }
}

setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12/phyper")
write.csv(results_table,file="Results_quartiles_phyper_pvalues.csv")