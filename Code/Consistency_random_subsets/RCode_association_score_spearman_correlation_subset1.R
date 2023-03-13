library(dplyr)

#Select cell types
cell_type_list=list(
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
  "Ast",
  "Oli",                       
  "OPC",
  "Mic"
)

for (c in 1:length(cell_type_list)) {
  cell_type=cell_type_list[[c]]


working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/gpath")
setwd(working_directory)

file_to_load=paste0(cell_type,"_gpath.csv")
topgenes=read.csv(file_to_load,header=TRUE,row.names=1)
Marker_genes_full=topgenes
Marker_genes_full_significant=subset(Marker_genes_full,p_adj.loc < 0.05)
Marker_genes_full_significant[,"score"]=(-log10(Marker_genes_full_significant[,"p_adj.loc"]))*((Marker_genes_full_significant[,"logFC"])/(abs(Marker_genes_full_significant[,"logFC"])))
Sign_genes_one=Marker_genes_full_significant[,c("gene","score")]
Sign_genes1=Marker_genes_full_significant[,c("gene"),drop=FALSE]

for (v in 1:length(cell_type_list)) {
  cell_type2=cell_type_list[[v]]
  
working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Revision_Dec2022/RUVr_random_subsets/Two_subsets/Results/gpath/random_subset_1")
setwd(working_directory)
file_to_load2=paste0(cell_type2,"_gpath.csv")
topgenes2=read.csv(file_to_load2,header=TRUE,row.names=1)
topgenes2[,"score"]=(-log10(topgenes2[,"p_adj.loc"]))*((topgenes2[,"logFC"])/(abs(topgenes2[,"logFC"])))
data_new_selected2=left_join(Sign_genes1,topgenes2,by="gene")
Sign_genes_two=data_new_selected2[,c("gene","score")]
names(Sign_genes_two)[names(Sign_genes_two) == "score"] <- "score2"
  
combined=left_join(Sign_genes_one,Sign_genes_two,by="gene")
  
  OK <- complete.cases(combined$score,combined$score2)
  x <- combined$score[OK]
  y <- combined$score2[OK]
  m <- length(x)
  if (m < 3L){
    results_p <- "NA"
    results_p = as.data.frame(results_p)
    rownames(results_p)=cell_type
    colnames(results_p)=cell_type2
    
    results_r <- "NA"
    results_r = as.data.frame(results_r)
    rownames(results_r)=cell_type
    colnames(results_r)=cell_type2
    
    if (v==1){
      results_df_p=results_p
      results_df_r=results_r
    }
    if (v>1){
      results_df_p=cbind(results_df_p,results_p)
      results_df_r=cbind(results_df_r,results_r)
    }
  } else {
    cor_results=cor.test(combined$score,combined$score2,alternative = c("two.sided"),method=c("spearman"))
    results_p=as.data.frame(-log10(cor_results[["p.value"]]))
    rownames(results_p)=cell_type
    colnames(results_p)=cell_type2
    
    results_r=as.data.frame(cor_results[["estimate"]])
    rownames(results_r)=cell_type
    colnames(results_r)=cell_type2
    
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
if (c==1){
  results_final_p=results_df_p
  results_final_r=results_df_r
}
if (c>1){
  results_final_p=rbind(results_final_p,results_df_p)
  results_final_r=rbind(results_final_r,results_df_r)
}
}
results_final_p[sapply(results_final_p, is.infinite)]  <- 300
setwd("D:/PFC_427_WS/Vanshika/muscat/Revision_Dec2022/Consistency_random_subsets//Results_effect_size")
output_file_r=paste0(cell_type,"_r.csv")
output_file_p=paste0(cell_type,"_p.csv")
write.csv(results_final_p,file="Association_scores_spearman_subset1_p.csv")
write.csv(results_final_r,file="Association_scores_spearman_subset1_r.csv")


