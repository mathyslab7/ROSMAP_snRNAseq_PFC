library(ggplot2)
library(Seurat)
library(readxl)  # install.packages("readxl") or install.packages("tidyverse")
library(plyr)
library(tibble)
library(dplyr)
library(data.table)
library(stringr)
library(qpcR)

setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex/SuppData")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")
corum_results2 <- subset(corum_results2, (ComplexID %in% c(165,166,5464,282,2839,3197,786,7313,7567,6529,6166,7459,7220,6994)))
data_split=split(corum_results2,corum_results2$ComplexID)
data_split2=lapply(data_split,"[", ,"subunits.Gene.name.")
data_split2=lapply(data_split2,function(z){z[z!=""]})

data_split3=list()
for(i in 1:length(data_split2)){
  data3=strsplit(data_split2[[i]],";", fixed = TRUE)
  data4=unlist(data3)
  data4=unique(data4)
  data_split3[[i]]=data4
  name=names(data_split2)[i]
  names(data_split3)[i]=name
}

data_split3=data_split3[lengths(data_split3) > 1]

Marker_genes_list=data_split3

###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "MT_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="MT_metadata.csv")

###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "AG_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="AG_metadata.csv")

###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "EC_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="EC_metadata.csv")


###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "HC_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="HC_metadata.csv")


###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "TH_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="TH_metadata.csv")


###################################################
###################################################

setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc <- readRDS(file = "PFC_merge_major_cell_type.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

all_genes=rownames(pbmc)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}

Marker_genes_list_pruned=Marker_genes_list_pruned[lengths(Marker_genes_list_pruned) > 1]

for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  pbmc <- AddModuleScore(object = pbmc, features = genes.for.scoring, name = name) 
}


setwd("G:/PFC_429_final_Spring2021/Vanshika/CORUM_complex_brainregion/Module_scores")
write.csv(pbmc@meta.data,file="PFC_metadata.csv")

###################################################
###################################################
