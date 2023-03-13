library(dplyr)
library(GeneOverlap)
data(GeneOverlap)


#Select variables
variable_list=list(
  "tangles"
)

for (k in 1:length(variable_list)) {
  variable=variable_list[[k]]
  
  working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",variable)
  setwd(working_directory)
  
  celltype = c("Exc L2-3 CBLN2 LINC02306","Exc L3-4 RORB CUX2","Exc L3-5 RORB PLCH1","Exc L4-5 RORB IL1RAPL2","Exc L4-5 RORB GABRG1",    
               "Exc L5-6 RORB LINC02196","Exc L6 THEMIS NFIA","Exc L5/6 IT Car3","Exc L5 ET","Exc L5/6 NP","Exc L6 CT","Exc L6b",
               "Inh L1-6 LAMP5 CA13","Inh LAMP5 NRG1 (Rosehip)","Inh PVALB SULF1","Inh PVALB CA8 (Chandelier)","Inh PVALB HTR4",
               "Inh CUX2 MSR1","Inh L3-5 SST MAFB")
  
  
  #List significantly upregulated genes
  
  Marker_genes_list<-list()
  for (i in 1:length(celltype)) {
    ct1=celltype[i]
    ct1=gsub("/","_",ct1)
    input_file=paste0(ct1,"_",variable,".csv")
    topgenes=read.csv(input_file,header=TRUE,row.names=1)
    Marker_genes_full=topgenes
    Marker_genes_full_significant=subset(Marker_genes_full,p_adj.loc < 0.05)
    Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC < 0)
    Marker_genes=Marker_genes_full_significant[,1:2]
    Marker_genes[] <- lapply(Marker_genes, as.character)
    if (nrow(Marker_genes)>0){
      Marker_genes_split=split(Marker_genes, Marker_genes$cluster_id)
      
      name=names(Marker_genes_split)[1]
      Marker_genes_list[[i]]<-Marker_genes_split[[1]]$gene
      names(Marker_genes_list)[i]<-name
    }
  }
  
  #Remove empty elements
  Marker_genes_list=Marker_genes_list[lengths(Marker_genes_list) != 0]
  
  for (c in 1:length(celltype)) {
    ct1=celltype[[c]]
    ct1_name=gsub("/","_",ct1)
    for (b in 1:length(celltype)){
      ct2=celltype[[b]]
      ct2_name=gsub("/","_",ct2)
      
      input_file1=paste0(ct1_name,"_",variable,".csv")
      input_file2=paste0(ct2_name,"_",variable,".csv")
      
      topgenes1=read.csv(input_file1,header=TRUE,row.names=1)
      
      topgenes2=read.csv(input_file2,header=TRUE,row.names=1)
      
      All_genes=full_join(topgenes1,topgenes2,by="gene")
      
      genome_size=nrow(All_genes)
      
      gom.obj <- newGOM(Marker_genes_list, Marker_genes_list,genome.size = genome_size)
      GeneOverlap_Matrix=getMatrix(gom.obj, name="pval")
      GeneOverlap_Matrix_df=as.data.frame(GeneOverlap_Matrix)
      
      myData=GeneOverlap_Matrix_df[ct1,ct2,drop=FALSE]
      
      if (b==1){
        results_int=myData
      }
      if(b>1){
        results_int=cbind(results_int,myData)
      }
    }
    if (c==1){
      results_table=results_int
    }
    if (c>1){
      results_table=rbind(results_table,results_int)
    }
  }
  
  n_rows=nrow(results_table)
  n_columns=ncol(results_table)
  factor=n_rows*n_columns
  pvalue_adjusted_matrix <- apply(results_table,2,p.adjust,method="bonferroni", n = factor)
  pvalue_adjusted_matrix_minusLog10=-log10(pvalue_adjusted_matrix)
  pvalue_adjusted_matrix_minusLog10[is.infinite(pvalue_adjusted_matrix_minusLog10)] <- 300
  
  setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap_v2/Overlap_celltypes/Results")
  output_file=paste0("MinusLog10BF_down_",variable,".csv")
  write.csv(pvalue_adjusted_matrix_minusLog10,file=output_file)
  
}
