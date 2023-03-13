library(Seurat)
library(dplyr)
library(Hmisc)
###########################################################
complexes=c("X1651","X1661","X54641")

#AG
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
metadata=read.csv("AG_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

metadata2=metadata[metadata$cell_type_high_resolution %in% c("Ast DPP10","Ast GRM3","Ast DCLK1","Exc L2-3 CBLN2 LINC02306","Exc L4-5 RORB IL1RAPL2"      
                                                 ,"Exc L3-4 RORB CUX2","Exc L4-5 RORB GABRG1","Exc L5/6 NP"                 
                                                 ,"Exc L6 THEMIS NFIA","Exc NRGN","Exc L6 CT"                   
                                                 ,"Exc L6b","Exc L3-5 RORB PLCH1","Exc L5-6 RORB LINC02196"     
                                                 ,"Exc L5/6 IT Car3","Exc L5 ET","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5"                  
                                                 ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                                                 ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                                                 ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                                                 ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                                                 ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                                                 ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                                                 ,"Inh LAMP5 RELN","Oli OPALIN","Oli RASGRF1","Per","End","Fib"),]

HM=metadata2

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

metadata3=HM3

metadata4=metadata3[metadata3$major_cell_type %in% c("Exc","Inh","Ast","Oli","OPC","Mic"),]
metadata5=metadata4[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "AG_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"

pbmc=subset(pbmc,idents = c("Ast DPP10","Ast GRM3","Ast DCLK1","Exc L2-3 CBLN2 LINC02306","Exc L4-5 RORB IL1RAPL2"      
                            ,"Exc L3-4 RORB CUX2","Exc L4-5 RORB GABRG1","Exc L5/6 NP"                 
                            ,"Exc L6 THEMIS NFIA","Exc NRGN","Exc L6 CT"                   
                            ,"Exc L6b","Exc L3-5 RORB PLCH1","Exc L5-6 RORB LINC02196"     
                            ,"Exc L5/6 IT Car3","Exc L5 ET","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5"                  
                            ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                            ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                            ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                            ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                            ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                            ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                            ,"Inh LAMP5 RELN","Oli OPALIN","Oli RASGRF1","Per","End","Fib"))

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
pbmc=subset(pbmc, idents=c("Exc","Inh","Ast","Oli","OPC","Mic"))

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion")
Union_genes=readRDS("Union_Exc_Oli_gpath_up.rds")
genes=Union_genes$gene

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata5, major_cell_type==celltype)
  for (j in 1:length(genes)){
    for (k in 1:length(complexes)){
    complex=complexes[[k]]
    gene=genes[[j]]
    expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
    expression_data[,"Cell_barcode"]=rownames(expression_data)
    combined=left_join(metadata_new,expression_data,by="Cell_barcode")
    results=cor.test(combined[,gene], combined[,complex],
             alternative = c("two.sided"),
             method = c("spearman"))
    results_rho=as.data.frame(results[["estimate"]][["rho"]])
    rownames(results_rho)=gene
    colnames(results_rho)=complex
    results_P=as.data.frame(results[["p.value"]])
    rownames(results_P)=gene
    colnames(results_P)=complex
    if (k==1){
      results_table_rho=results_rho
      results_table_P=results_P
    }
    if (k>1){
      results_table_rho=cbind(results_table_rho,results_rho)
      results_table_P=cbind(results_table_P,results_P)
    }
    }
    if (j==1){
      results_table_rho_final=results_table_rho
      results_table_P_final=results_table_P
    }
    if (j>1){
      results_table_rho_final=rbind(results_table_rho_final,results_table_rho)
      results_table_P_final=rbind( results_table_P_final,results_table_P)
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0("AG_",celltype,"_results_table_rho.csv")
  filename_P=paste0("AG_",celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Correlation_results")
  write.csv(results_table_rho_final,file=filename_r)
  write.csv(results_table_P_final,file=filename_P)
}

###########################################################


#MT
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
metadata=read.csv("MT_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

metadata2=metadata[metadata$cell_type_high_resolution %in% c("Ast DPP10","Ast GRM3","Exc L2-3 CBLN2 LINC02306","Exc L4-5 RORB IL1RAPL2"      
                                                             ,"Exc L3-4 RORB CUX2","Exc L4-5 RORB GABRG1","Exc L5/6 NP"                 
                                                             ,"Exc L6 THEMIS NFIA","Exc NRGN","Exc L6 CT"                   
                                                             ,"Exc L6b","Exc L3-5 RORB PLCH1","Exc L5-6 RORB LINC02196"     
                                                             ,"Exc L5/6 IT Car3","Exc L5 ET","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                                                             ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                                                             ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                                                             ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                                                             ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                                                             ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                                                             ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                                                             ,"Inh LAMP5 RELN","Inh L6 SST NPY","Oli OPALIN","Oli RASGRF1","Per","End","Fib"),]

HM=metadata2

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

metadata3=HM3

metadata4=metadata3[metadata3$major_cell_type %in% c("Exc","Inh","Ast","Oli","OPC","Mic"),]
metadata5=metadata4[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "MT_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"

pbmc=subset(pbmc,idents = c("Ast DPP10","Ast GRM3","Exc L2-3 CBLN2 LINC02306","Exc L4-5 RORB IL1RAPL2"      
                            ,"Exc L3-4 RORB CUX2","Exc L4-5 RORB GABRG1","Exc L5/6 NP"                 
                            ,"Exc L6 THEMIS NFIA","Exc NRGN","Exc L6 CT"                   
                            ,"Exc L6b","Exc L3-5 RORB PLCH1","Exc L5-6 RORB LINC02196"     
                            ,"Exc L5/6 IT Car3","Exc L5 ET","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                            ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                            ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                            ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                            ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                            ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                            ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                            ,"Inh LAMP5 RELN","Inh L6 SST NPY","Oli OPALIN","Oli RASGRF1","Per","End","Fib"))

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
pbmc=subset(pbmc, idents=c("Exc","Inh","Ast","Oli","OPC","Mic"))

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion")
Union_genes=readRDS("Union_Exc_Oli_gpath_up.rds")
genes=Union_genes$gene

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata5, major_cell_type==celltype)
  for (j in 1:length(genes)){
    for (k in 1:length(complexes)){
      complex=complexes[[k]]
      gene=genes[[j]]
      expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
      expression_data[,"Cell_barcode"]=rownames(expression_data)
      combined=left_join(metadata_new,expression_data,by="Cell_barcode")
      results=cor.test(combined[,gene], combined[,complex],
                       alternative = c("two.sided"),
                       method = c("spearman"))
      results_rho=as.data.frame(results[["estimate"]][["rho"]])
      rownames(results_rho)=gene
      colnames(results_rho)=complex
      results_P=as.data.frame(results[["p.value"]])
      rownames(results_P)=gene
      colnames(results_P)=complex
      if (k==1){
        results_table_rho=results_rho
        results_table_P=results_P
      }
      if (k>1){
        results_table_rho=cbind(results_table_rho,results_rho)
        results_table_P=cbind(results_table_P,results_P)
      }
    }
    if (j==1){
      results_table_rho_final=results_table_rho
      results_table_P_final=results_table_P
    }
    if (j>1){
      results_table_rho_final=rbind(results_table_rho_final,results_table_rho)
      results_table_P_final=rbind( results_table_P_final,results_table_P)
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0("MT_",celltype,"_results_table_rho.csv")
  filename_P=paste0("MT_",celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Correlation_results")
  write.csv(results_table_rho_final,file=filename_r)
  write.csv(results_table_P_final,file=filename_P)
}

###########################################################


#EC
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
metadata=read.csv("EC_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

metadata2=metadata[metadata$cell_type_high_resolution %in% c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc RELN COL5A2","Exc DLC1 SNTG2",
                                                             "Exc TOX3 TTC6","Exc RELN GPC5","Exc TOX3 INO80D","Exc AGBL1 GPC5","Exc SOX11 NCKAP5",
                                                             "Exc TOX3 POSTN","Exc COL25A1 SEMA3D","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                                                             ,"Inh GPC5 RIT2","Inh LAMP5 NRG1 (Rosehip)","Inh SGCD PDE3A",              
                                                             "Inh VIP CLSTN2","Inh RYR3 TSHZ2","Inh SORCS1 TTN",             
                                                             "Inh L1-6 LAMP5 CA13","Inh CUX2 MSR1","Inh VIP ABI3BP",              
                                                             "Inh FBN2 EPB41L4A","Inh PVALB HTR4","Inh L3-5 SST MAFB",           
                                                             "Inh PAX6 RELN","Inh PVALB SULF1","Inh VIP TSHZ2",               
                                                             "Inh LAMP5 RELN","Inh ALCAM TRPM3","Inh PVALB CA8 (Chandelier)",  
                                                             "Inh ENOX2 SPHKAP","Inh L6 SST NPY","Inh VIP THSD7B",             
                                                             "Inh PTPRK FAM19A1","Oli OPALIN","Oli RASGRF1","Per","End","Fib","SMC"),]

HM=metadata2

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

metadata3=HM3

metadata4=metadata3[metadata3$major_cell_type %in% c("Exc","Inh","Ast","Oli","OPC","Mic"),]
metadata5=metadata4[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "EC_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"

pbmc=subset(pbmc,idents = c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc RELN COL5A2","Exc DLC1 SNTG2",
                            "Exc TOX3 TTC6","Exc RELN GPC5","Exc TOX3 INO80D","Exc AGBL1 GPC5","Exc SOX11 NCKAP5",
                            "Exc TOX3 POSTN","Exc COL25A1 SEMA3D","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                            ,"Inh GPC5 RIT2","Inh LAMP5 NRG1 (Rosehip)","Inh SGCD PDE3A",              
                            "Inh VIP CLSTN2","Inh RYR3 TSHZ2","Inh SORCS1 TTN",             
                            "Inh L1-6 LAMP5 CA13","Inh CUX2 MSR1","Inh VIP ABI3BP",              
                            "Inh FBN2 EPB41L4A","Inh PVALB HTR4","Inh L3-5 SST MAFB",           
                            "Inh PAX6 RELN","Inh PVALB SULF1","Inh VIP TSHZ2",               
                            "Inh LAMP5 RELN","Inh ALCAM TRPM3","Inh PVALB CA8 (Chandelier)",  
                            "Inh ENOX2 SPHKAP","Inh L6 SST NPY","Inh VIP THSD7B",             
                            "Inh PTPRK FAM19A1","Oli OPALIN","Oli RASGRF1","Per","End","Fib","SMC"))

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
pbmc=subset(pbmc, idents=c("Exc","Inh","Ast","Oli","OPC","Mic"))

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion")
Union_genes=readRDS("Union_Exc_Oli_gpath_up.rds")
genes=Union_genes$gene

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata5, major_cell_type==celltype)
  for (j in 1:length(genes)){
    for (k in 1:length(complexes)){
      complex=complexes[[k]]
      gene=genes[[j]]
      expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
      expression_data[,"Cell_barcode"]=rownames(expression_data)
      combined=left_join(metadata_new,expression_data,by="Cell_barcode")
      results=cor.test(combined[,gene], combined[,complex],
                       alternative = c("two.sided"),
                       method = c("spearman"))
      results_rho=as.data.frame(results[["estimate"]][["rho"]])
      rownames(results_rho)=gene
      colnames(results_rho)=complex
      results_P=as.data.frame(results[["p.value"]])
      rownames(results_P)=gene
      colnames(results_P)=complex
      if (k==1){
        results_table_rho=results_rho
        results_table_P=results_P
      }
      if (k>1){
        results_table_rho=cbind(results_table_rho,results_rho)
        results_table_P=cbind(results_table_P,results_P)
      }
    }
    if (j==1){
      results_table_rho_final=results_table_rho
      results_table_P_final=results_table_P
    }
    if (j>1){
      results_table_rho_final=rbind(results_table_rho_final,results_table_rho)
      results_table_P_final=rbind( results_table_P_final,results_table_P)
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0("EC_",celltype,"_results_table_rho.csv")
  filename_P=paste0("EC_",celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Correlation_results")
  write.csv(results_table_rho_final,file=filename_r)
  write.csv(results_table_P_final,file=filename_P)
}

###########################################################


#HC
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
metadata=read.csv("HC_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

metadata2=metadata[metadata$cell_type_high_resolution %in% c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc TRPC6 ANO2","Exc GRIK1 CTXND1 (Subiculum)","Exc ZNF385D COL24A1","Exc COBLL1 UST"
                                                             ,"CA2, CA3 pyramidal cells","DG granule cells","CA1 pyramidal cells" ,"Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                                                             ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                                                             ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                                                             ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                                                             ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                                                             ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                                                             ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                                                             ,"Inh LAMP5 RELN","Inh L6 SST NPY","Oli OPALIN","Oli RASGRF1","Per","End","Fib","SMC","Epd","CPEC"),]

HM=metadata2

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

metadata3=HM3

metadata4=metadata3[metadata3$major_cell_type %in% c("Exc","Inh","Ast","Oli","OPC","Mic"),]
metadata5=metadata4[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "HC_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"

pbmc=subset(pbmc,idents = c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc TRPC6 ANO2","Exc GRIK1 CTXND1 (Subiculum)","Exc ZNF385D COL24A1","Exc COBLL1 UST"
                            ,"CA2, CA3 pyramidal cells","DG granule cells","CA1 pyramidal cells" ,"Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                            ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                            ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                            ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                            ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                            ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                            ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                            ,"Inh LAMP5 RELN","Inh L6 SST NPY","Oli OPALIN","Oli RASGRF1","Per","End","Fib","SMC","Epd","CPEC"))

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
pbmc=subset(pbmc, idents=c("Exc","Inh","Ast","Oli","OPC","Mic"))

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion")
Union_genes=readRDS("Union_Exc_Oli_gpath_up.rds")
genes=Union_genes$gene

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata5, major_cell_type==celltype)
  for (j in 1:length(genes)){
    for (k in 1:length(complexes)){
      complex=complexes[[k]]
      gene=genes[[j]]
      expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
      expression_data[,"Cell_barcode"]=rownames(expression_data)
      combined=left_join(metadata_new,expression_data,by="Cell_barcode")
      results=cor.test(combined[,gene], combined[,complex],
                       alternative = c("two.sided"),
                       method = c("spearman"))
      results_rho=as.data.frame(results[["estimate"]][["rho"]])
      rownames(results_rho)=gene
      colnames(results_rho)=complex
      results_P=as.data.frame(results[["p.value"]])
      rownames(results_P)=gene
      colnames(results_P)=complex
      if (k==1){
        results_table_rho=results_rho
        results_table_P=results_P
      }
      if (k>1){
        results_table_rho=cbind(results_table_rho,results_rho)
        results_table_P=cbind(results_table_P,results_P)
      }
    }
    if (j==1){
      results_table_rho_final=results_table_rho
      results_table_P_final=results_table_P
    }
    if (j>1){
      results_table_rho_final=rbind(results_table_rho_final,results_table_rho)
      results_table_P_final=rbind( results_table_P_final,results_table_P)
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0("HC_",celltype,"_results_table_rho.csv")
  filename_P=paste0("HC_",celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Correlation_results")
  write.csv(results_table_rho_final,file=filename_r)
  write.csv(results_table_P_final,file=filename_P)
}

###########################################################


#TH
###########################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
metadata=read.csv("TH_metadata.csv",row.names = 1)
metadata[,"Cell_barcode"]=rownames(metadata)

metadata2=metadata[metadata$cell_type_high_resolution %in% c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc NXPH1 RNF220","Exc VAT1L ERBB4","Exc NRGN",             
                                                             "Mic P2RY12","Mic TPT1","Mic DUSP1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                                                             ,"Inh LAMP5 NRG1 (Rosehip)","Inh MEIS2 FOXP2","Oli OPALIN","Oli RASGRF1","Per","End","Fib","Epd"),]

HM=metadata2

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

metadata3=HM3

metadata4=metadata3[metadata3$major_cell_type %in% c("Exc","Inh","Ast","Oli","OPC","Mic"),]
metadata5=metadata4[,c("Cell_barcode","major_cell_type","X1651","X1661","X54641")]

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "TH_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"

pbmc=subset(pbmc,idents = c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc NXPH1 RNF220","Exc VAT1L ERBB4","Exc NRGN",             
                            "Mic P2RY12","Mic TPT1","Mic DUSP1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                            ,"Inh LAMP5 NRG1 (Rosehip)","Inh MEIS2 FOXP2","Oli OPALIN","Oli RASGRF1","Per","End","Fib","Epd"))

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
pbmc=subset(pbmc, idents=c("Exc","Inh","Ast","Oli","OPC","Mic"))

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion")
Union_genes=readRDS("Union_Exc_Oli_gpath_up.rds")
genes=Union_genes$gene

object_list=SplitObject(pbmc, split.by = "major_cell_type")
rm(pbmc)

for (i in 1:length(object_list)){
  celltype=names(object_list)[i]
  metadata_new=subset(metadata5, major_cell_type==celltype)
  for (j in 1:length(genes)){
    for (k in 1:length(complexes)){
      complex=complexes[[k]]
      gene=genes[[j]]
      expression_data=FetchData(object = object_list[[i]], vars=gene, slot = "data")
      expression_data[,"Cell_barcode"]=rownames(expression_data)
      combined=left_join(metadata_new,expression_data,by="Cell_barcode")
      results=cor.test(combined[,gene], combined[,complex],
                       alternative = c("two.sided"),
                       method = c("spearman"))
      results_rho=as.data.frame(results[["estimate"]][["rho"]])
      rownames(results_rho)=gene
      colnames(results_rho)=complex
      results_P=as.data.frame(results[["p.value"]])
      rownames(results_P)=gene
      colnames(results_P)=complex
      if (k==1){
        results_table_rho=results_rho
        results_table_P=results_P
      }
      if (k>1){
        results_table_rho=cbind(results_table_rho,results_rho)
        results_table_P=cbind(results_table_P,results_P)
      }
    }
    if (j==1){
      results_table_rho_final=results_table_rho
      results_table_P_final=results_table_P
    }
    if (j>1){
      results_table_rho_final=rbind(results_table_rho_final,results_table_rho)
      results_table_P_final=rbind( results_table_P_final,results_table_P)
    }
  }
  celltype=gsub("/","_",celltype)
  filename_r=paste0("TH_",celltype,"_results_table_rho.csv")
  filename_P=paste0("TH_",celltype,"_results_table_P.csv")
  setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Correlation_results")
  write.csv(results_table_rho_final,file=filename_r)
  write.csv(results_table_P_final,file=filename_P)
}

###########################################################