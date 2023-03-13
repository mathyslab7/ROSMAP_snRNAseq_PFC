library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_40_percent.rds")

DefaultAssay(pbmc)="RNA"

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Add major cell type info
Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic",cell_type_high_resolution)~"Mic",
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

#Subset quartiles to be compared
Idents(pbmc)<-'major_cell_type'
pbmc=subset(pbmc,idents=c('Exc'))

average_expression=AverageExpression(pbmc)

setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Average_expression")
write.csv(average_expression[["RNA"]],file="Average_expression_Exc_neurons.csv")