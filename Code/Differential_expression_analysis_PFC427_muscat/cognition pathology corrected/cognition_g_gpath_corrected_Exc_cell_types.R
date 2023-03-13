


library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)


##########################################################################################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set1_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "cell_type_high_resolution", # subpopulation assignments
                gid = "AD",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

rm(pbmc)
rm(sce)
rm(pbmc.sce)

##############################################
formula = ~cogn_global_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc1_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc1_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set2_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "cell_type_high_resolution", # subpopulation assignments
                gid = "AD",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

rm(pbmc)
rm(sce)
rm(pbmc.sce)


##############################################
formula = ~cogn_global_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc2_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc2_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set3_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "cell_type_high_resolution", # subpopulation assignments
                gid = "AD",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc3_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  name=gsub("/","_",name)
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/Exc3_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_gpath_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  name=gsub("/","_",name)
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

