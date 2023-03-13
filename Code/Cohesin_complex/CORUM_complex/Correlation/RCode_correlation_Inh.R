library(Seurat)
library(dplyr)
library(Hmisc)
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex/Module_scores")
metadata=read.csv("In_integrated_batch_3000_module_scores_predicted_final_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)
metadata=metadata[,c("Cell_barcode","cell_type_high_resolution","X1651","X1661","X54641")]

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "In_integrated_batch_3000_module_scores_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
genes=rownames(pbmc)

object_list=SplitObject(pbmc, split.by = "cell_type_high_resolution")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata, cell_type_high_resolution==celltype)
  for (j in 1:length(genes)){
    gene=genes[[j]]
    expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
    expression_data[,"Cell_barcode"]=rownames(expression_data)
    combined=left_join(metadata_new,expression_data,by="Cell_barcode")
    results=rcorr(as.matrix(combined[,3:6]), type = c("spearman"))
    results_r=as.data.frame(results[["r"]])
    results_P=as.data.frame(results[["P"]])
    if (j==1){
      results_table_r=results_r[1:3,gene,drop=FALSE]
      results_table_P=results_P[1:3,gene,drop=FALSE]
    }
    if (j>1){
      results_table_r=cbind(results_table_r,results_r[1:3,gene,drop=FALSE])
      results_table_P=cbind(results_table_P,results_P[1:3,gene,drop=FALSE])
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0(celltype,"_results_table_r.csv")
  filename_P=paste0(celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex/Correlation_results")
  write.csv(results_table_r,file=filename_r)
  write.csv(results_table_P,file=filename_P)
}

###########################################################