library(dplyr)

#Select cell types
cell_type_list=list(
  "Ast CHI3L1",                
  "Ast DPP10",                 
  "Ast GRM3",
  "Ast",
  "CAMs",                     
  "End",                       
  "Exc L2-3 CBLN2 LINC02306",  
  "Exc L3-4 RORB CUX2",        
  "Exc L3-5 RORB PLCH1",       
  "Exc L4-5 RORB GABRG1",      
  "Exc L4-5 RORB IL1RAPL2",    
  "Exc L5-6 RORB LINC02196",   
  "Exc L5 ET",                 
  "Exc L5_6 IT Car3",          
  "Exc L5_6 NP",               
  "Exc L6 CT",                 
  "Exc L6 THEMIS NFIA",        
  "Exc L6b",                   
  "Exc NRGN",                  
  "Exc RELN CHD7",             
  "Fib FLRT2",                 
  #"Fib SLC4A4",                
  "Inh ALCAM TRPM3",           
  "Inh CUX2 MSR1",             
  "Inh ENOX2 SPHKAP",          
  "Inh FBN2 EPB41L4A",         
  "Inh GPC5 RIT2",             
  #"Inh L1-2 PAX6 SCGN",        
  "Inh L1-6 LAMP5 CA13",       
  "Inh L1 PAX6 CA4",           
  "Inh L3-5 SST MAFB",         
  "Inh L5-6 PVALB STON2",      
  "Inh L5-6 SST TH",           
  #"Inh L6 SST NPY",            
  #"Inh LAMP5 NRG1 (Rosehip)",  
  "Inh LAMP5 RELN",            
  "Inh PTPRK FAM19A1",         
  #"Inh PVALB CA8 (Chandelier)",
  "Inh PVALB HTR4",            
  "Inh PVALB SULF1",           
  "Inh RYR3 TSHZ2",            
  "Inh SGCD PDE3A",            
  "Inh SORCS1 TTN",            
  "Inh VIP ABI3BP",            
  "Inh VIP CLSTN2",            
  "Inh VIP THSD7B",            
  "Inh VIP TSHZ2",             
  #"Mic MKI67",                 
  "Mic P2RY12",                
  "Mic TPT1",                  
  "Mic",                       
  "Oli",                       
  "OPC",                       
  "Per"                   
  #"SMC",                       
  #"T cells"
)

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

for (c in 1:length(cell_type_list)) {
  cell_type=cell_type_list[[c]]

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

#cell_type="Exc L2-3 CBLN2 LINC02306"
myFiles=character()
for (k in 1:length(variable_list)) {
variable=variable_list[[k]]

folder_variable=variable_list[[k]]

working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",folder_variable)
setwd(working_directory)

pattern=paste0(cell_type,"_",variable,".csv")
#pattern=paste0(cell_type,"_*")
#Select major cell type
#Reading in all the ".csv" files from the working directory starting with "Exc"
myFiles1 <- list.files(pattern=pattern)
#myFiles2 <- list.files(pattern="Inh*")
myFiles <-append(myFiles, myFiles1)
}

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
Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC > 0)
}
if(any(grepl(variable,Group_2))){
  Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC < 0)
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


setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap/Overlap_variables/path_upreg_cog_downreg")
file_name=paste0("GeneOverlap_BFminusLOG10_",cell_type,".csv")
write.csv(logp_df, file=file_name)
}
