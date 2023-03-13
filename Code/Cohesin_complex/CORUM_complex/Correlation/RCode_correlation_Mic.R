library(Seurat)
library(dplyr)
library(Hmisc)
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex/Module_scores")
metadata=read.csv("Immune_integrated_batch_5000_predicted_final_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

HM=metadata

HM3 = HM %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh ",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast ",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli ",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC ",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic ",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))

metadata=HM3
metadata=metadata[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Immune_integrated_batch_5000_predicted_final.rds")

Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh ",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast ",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli ",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC ",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic ",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))
HM4=HM3[,-1]
HM5=HM4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,HM5)
Idents(pbmc) <- 'major_cell_type'
pbmc=subset(pbmc, idents="Mic")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
genes=rownames(pbmc)

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata, major_cell_type==celltype)
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