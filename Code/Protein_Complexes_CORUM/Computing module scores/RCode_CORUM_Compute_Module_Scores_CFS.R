library(ggplot2)
library(Seurat)
library(dplyr)

###################################

#Ast
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Ast_integrated_batch_5000_leiden_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Ast_integrated_batch_5000_leiden_predicted_final_metadata.csv")

###################################

#Exc1
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Exc_raw_set1_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Exc_raw_set1_final_metadata.csv")

###################################

#Exc2
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Exc_raw_set2_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Exc_raw_set2_final_metadata.csv")

###################################

#Exc3
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Exc_raw_set3_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Exc_raw_set3_final_metadata.csv")

###################################

#Inh
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "In_integrated_batch_3000_module_scores_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="In_integrated_batch_3000_module_scores_predicted_final_metadata.csv")

###################################

#Oli
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Oli_integrated_batch3000_leiden_res02_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Oli_integrated_batch3000_leiden_res02_final_metadata.csv")

###################################

#OPC
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Opc_integrated_batch_5000_leiden_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Opc_integrated_batch_5000_leiden_predicted_final_metadata.csv")

###################################

#Immune
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Immune_integrated_batch_5000_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Immune_integrated_batch_5000_predicted_final_metadata.csv")

###################################

#Vasc
setwd("F:/PFC_429_final_Spring2021/Data")
pbmc <- readRDS(file = "Vasc_integrated_batch_5000_predicted_final.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Supp_data")
corum_results=read.csv("coreComplexes.csv")
corum_results2=subset(corum_results, Organism=="Human")

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

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
write.csv(pbmc@meta.data,file="Vasc_integrated_batch_5000_predicted_final_metadata.csv")

###################################