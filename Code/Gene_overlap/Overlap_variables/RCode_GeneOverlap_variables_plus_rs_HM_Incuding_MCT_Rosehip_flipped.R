library(dplyr)

#Select cell types
#"Inh LAMP5 NRG1 (Rosehip)",
#"Inh PVALB CA8 (Chandelier)",


# Select variables' group (directionality)
Group_1=list(
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
  "bradysc_lv",
  "gaitsc_lv",
  "parksc_lv",
  "arteriol_scler",
  "caa_4gp",
  "cvda_4gp2",
  "ci_num2_gct",
  "ci_num2_gtt",
  "ci_num2_mct",
  "ci_num2_mtt",
  "diabetes_sr_rx_bl",
  "cancer_bl",
  "stroke_bl",
  "headinjrloc_bl",
  "heart_bl",
  "hypertension_bl",
  "tdp_st4",
  "dlbdx",
  "cogdx_stroke",
  "dxpark",
  "msex",
  "cogdx")

Group_2=list(
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
  "cognwo_random_slope"
)

#for (c in 1:length(cell_type_list)) {
#  cell_type=cell_type_list[[c]]

# Select variables
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
  "bradysc_lv",
  "gaitsc_lv",
  "parksc_lv",
  "arteriol_scler",
  "caa_4gp",
  "cvda_4gp2",
  "ci_num2_gct",
  "ci_num2_gtt",
  "ci_num2_mct",
  "ci_num2_mtt",
  "diabetes_sr_rx_bl",
  "cancer_bl",
  "stroke_bl",
  "headinjrloc_bl",
  "heart_bl",
  "hypertension_bl",
  "tdp_st4",
  "dlbdx",
  "cogdx_stroke",
  "dxpark",
  "msex",
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
  "cognwo_random_slope"
)





myFiles=c("Inh LAMP5 NRG1 (Rosehip)_amyloid.csv",
          "Inh LAMP5 NRG1 (Rosehip)_gpath.csv",
          "Inh LAMP5 NRG1 (Rosehip)_gpath_3neocort.csv",
          "Inh LAMP5 NRG1 (Rosehip)_nft.csv",
          "Inh LAMP5 NRG1 (Rosehip)_nft_mf.csv",
          "Inh LAMP5 NRG1 (Rosehip)_tangles.csv",
          "Inh LAMP5 NRG1 (Rosehip)_plaq_n.csv",
          "Inh LAMP5 NRG1 (Rosehip)_plaq_n_mf.csv",
          "Inh LAMP5 NRG1 (Rosehip)_plaq_d.csv",
          "Inh LAMP5 NRG1 (Rosehip)_plaq_d_mf.csv",
          "Inh LAMP5 NRG1 (Rosehip)_Apoe_e4.csv",
          "Inh LAMP5 NRG1 (Rosehip)_bradysc_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_gaitsc_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_parksc_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_arteriol_scler.csv",
          "Inh LAMP5 NRG1 (Rosehip)_caa_4gp.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cvda_4gp2.csv",
          "Inh LAMP5 NRG1 (Rosehip)_ci_num2_gct.csv",
          "Inh LAMP5 NRG1 (Rosehip)_ci_num2_gtt.csv",
          "Inh LAMP5 NRG1 (Rosehip)_ci_num2_mct.csv",
          "Inh LAMP5 NRG1 (Rosehip)_ci_num2_mtt.csv",
          "Inh LAMP5 NRG1 (Rosehip)_diabetes_sr_rx_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cancer_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_stroke_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_headinjrloc_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_heart_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_hypertension_bl.csv",
          "Inh LAMP5 NRG1 (Rosehip)_tdp_st4.csv",
          "Inh LAMP5 NRG1 (Rosehip)_dlbdx.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogdx_stroke.csv",
          "Inh LAMP5 NRG1 (Rosehip)_dxpark.csv",
          "Inh LAMP5 NRG1 (Rosehip)_msex.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogdx.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognition.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_global_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_ep_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_po_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_ps_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_se_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogn_wo_lv.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cogng_random_slope.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognep_random_slope.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognpo_random_slope.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognps_random_slope.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognse_random_slope.csv",
          "Inh LAMP5 NRG1 (Rosehip)_cognwo_random_slope.csv")


#List significantly upregulated genes
Marker_genes_list <- list()
for (i in 1:length(myFiles)) {
  variable=variable_list[[i]]
  folder_variable=variable_list[[i]]
  working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable)
  setwd(working_directory)
topgenes=read.csv(myFiles[i],header=TRUE,row.names=1)
Marker_genes_full=topgenes
Marker_genes_full_significant=subset(Marker_genes_full,p_adj.loc < 0.05)

if(any(grepl(variable,Group_1))){
Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC < 0)
}
if(any(grepl(variable,Group_2))){
  Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC > 0)
}

Marker_genes=Marker_genes_full_significant[,1:2]
Marker_genes[] <- lapply(Marker_genes, as.character)
if (nrow(Marker_genes)>0){
  Marker_genes_split=split(Marker_genes, Marker_genes$cluster_id)

  name=variable
  Marker_genes_list[[i]]<-Marker_genes_split[[1]]$gene
  names(Marker_genes_list)[i]<-name
  }
}

#Remove empty elements
Marker_genes_list=Marker_genes_list[lengths(Marker_genes_list) != 0]

variable=variable_list[[1]]
folder_variable=variable_list[[1]]
working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable)
setwd(working_directory)

All_genes<-data.frame()
topgenes=read.csv(myFiles[1],header=TRUE,row.names=1)
All_genes=topgenes

#for (j in 2:length(myFiles)) {
#  variable=variable_list[[j]]
#  working_directory=paste0("F:/PFC_429_final_Spring2021/Vanshika/Muscat_differential_gene_expression_analysis/",variable)
#  setwd(working_directory)
#  topgenes=read.csv(myFiles[j],header=TRUE,row.names=1)
#  All_genes=full_join(All_genes,topgenes,by="gene")
#}
genome_size=nrow(All_genes)

#setwd("F:/PFC_429_final_Spring2021/Vanshika/Muscat_differential_gene_expression_analysis/Gene_overlap/Overlap_Variables/variables_plus_rs/upreg_downreg")
library(GeneOverlap)
data(GeneOverlap)
gom.obj <- newGOM(Marker_genes_list, Marker_genes_list,genome.size = genome_size)
GeneOverlap_Matrix=getMatrix(gom.obj, name="pval")
GeneOverlap_Matrix_df=as.data.frame(GeneOverlap_Matrix)
GeneOverlap_Matrix_df[GeneOverlap_Matrix_df==0]<-NA
min=min(GeneOverlap_Matrix_df,na.rm=TRUE)

GeneOverlap_Matrix_df2=as.data.frame(GeneOverlap_Matrix)
GeneOverlap_Matrix_df2[GeneOverlap_Matrix_df2==0]<-min

#Bonferroni correction
factor=nrow(GeneOverlap_Matrix_df2)*ncol(GeneOverlap_Matrix_df2)
GeneOverlap_Matrix_df2=GeneOverlap_Matrix_df2*factor

#minus log10(pvalues)
logp_df=-log10(GeneOverlap_Matrix_df2)
logp_df[logp_df<0]<-0


setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap/Overlap_variables/path_downreg_cog_upreg")
file_name=paste0("GeneOverlap_BFminusLOG10_","Inh LAMP5 NRG1 (Rosehip)",".csv")
write.csv(logp_df, file=file_name)
#}
