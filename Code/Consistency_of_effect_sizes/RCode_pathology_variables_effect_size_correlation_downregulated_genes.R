library(dplyr)

#Select cell types
cell_type_list=list(
  "Ast",
  "Exc L2-3 CBLN2 LINC02306",  
  "Exc L3-4 RORB CUX2",        
  "Exc L3-5 RORB PLCH1",       
  "Exc L4-5 RORB GABRG1",      
  "Exc L4-5 RORB IL1RAPL2",    
  "Exc L5-6 RORB LINC02196",   
  "Exc L5_6 IT Car3",          
  "Exc L5_6 NP",               
  "Exc L6 CT",                 
  "Exc L6 THEMIS NFIA",        
  "Exc L6b",
  "Mic",                       
  "Oli",                       
  "OPC"
)

for (c in 1:length(cell_type_list)) {
  cell_type=cell_type_list[[c]]

# Select variables
variable_list=list(
  "amyloid",
  "gpath",
  "nft",
  "tangles",
  "plaq_n",
  "plaq_d")

for (k in 1:length(variable_list)) {
variable=variable_list[[k]]

folder_variable=variable_list[[k]]

working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable)
setwd(working_directory)

file_to_load=paste0(cell_type,"_",variable,".csv")
topgenes=read.csv(file_to_load,header=TRUE,row.names=1)
Marker_genes_full=topgenes
Marker_genes_full_significant=subset(Marker_genes_full,p_adj.loc < 0.05)
Marker_genes_full_significant=subset(Marker_genes_full,logFC < 0)
Sign_genes1=Marker_genes_full_significant[,c("gene"),drop=FALSE]

for (v in 1:length(variable_list)) {
  variable2=variable_list[[v]]
  folder_variable2=variable_list[[v]]
  working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable2)
  setwd(working_directory)
  file_to_load2=paste0(cell_type,"_",variable2,".csv")
  topgenes2=read.csv(file_to_load2,header=TRUE,row.names=1)
  Marker_genes_full=topgenes2
  Marker_genes_full_significant=subset(Marker_genes_full,p_adj.loc < 0.05)
  Marker_genes_full_significant=subset(Marker_genes_full,logFC < 0)
  Sign_genes2=Marker_genes_full_significant[,c("gene"),drop=FALSE]
  
  combined_list=full_join(Sign_genes1,Sign_genes2,by="gene")
  
  #load the data again
  working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable)
  setwd(working_directory)
  
  data_new=read.csv(file_to_load,header=TRUE,row.names=1)
  data_new_selected=left_join(combined_list,data_new,by="gene")
  Sign_genes_one=data_new_selected[,c("gene","logFC")]
  
  working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable2)
  setwd(working_directory)
  data_new2=read.csv(file_to_load2,header=TRUE,row.names=1)
  data_new_selected2=left_join(combined_list,data_new2,by="gene")
  Sign_genes_two=data_new_selected2[,c("gene","logFC")]
  names(Sign_genes_two)[names(Sign_genes_two) == "logFC"] <- "logFC2"
  
  combined=left_join(Sign_genes_one,Sign_genes_two,by="gene")
  
  OK <- complete.cases(combined$logFC,combined$logFC2)
  x <- combined$logFC[OK]
  y <- combined$logFC2[OK]
  m <- length(x)
  if (m < 3L){
    results_p <- "NA"
    results_p = as.data.frame(results_p)
    rownames(results_p)=variable
    colnames(results_p)=variable2
    
    results_r <- "NA"
    results_r = as.data.frame(results_r)
    rownames(results_r)=variable
    colnames(results_r)=variable2
    
    if (v==1){
      results_df_p=results_p
      results_df_r=results_r
    }
    if (v>1){
      results_df_p=cbind(results_df_p,results_p)
      results_df_r=cbind(results_df_r,results_r)
    }
  } else {
    cor_results=cor.test(combined$logFC,combined$logFC2,alternative = c("two.sided"),method=c("pearson"))
    results_p=as.data.frame(-log10(cor_results[["p.value"]]))
    rownames(results_p)=variable
    colnames(results_p)=variable2
    
    results_r=as.data.frame(cor_results[["estimate"]])
    rownames(results_r)=variable
    colnames(results_r)=variable2
    
    if (v==1){
      results_df_p=results_p
      results_df_r=results_r
    }
    if (v>1){
      results_df_p=cbind(results_df_p,results_p)
      results_df_r=cbind(results_df_r,results_r)
    }
  }
}
if (k==1){
  results_final_p=results_df_p
  results_final_r=results_df_r
}
if (k>1){
  results_final_p=rbind(results_final_p,results_df_p)
  results_final_r=rbind(results_final_r,results_df_r)
}
}
setwd("D:/PFC_427_WS/Vanshika/muscat/Revision_Dec2022/Consistency_of_effect_sizes/Results_downregulated_genes")
output_file_r=paste0(cell_type,"_r.csv")
output_file_p=paste0(cell_type,"_p.csv")
write.csv(results_final_p,file=output_file_p)
write.csv(results_final_r,file=output_file_r)
}

