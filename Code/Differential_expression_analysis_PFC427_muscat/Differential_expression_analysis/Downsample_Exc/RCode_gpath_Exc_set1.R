library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)


setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set1_final.rds")

#Split into individual cell types
Idents(pbmc)<-"cell_type_high_resolution"
pbmc.list <- SplitObject(pbmc, split.by = "cell_type_high_resolution")

cell_number_list=c(10,20,30,40,50,60,70,80,90,100)

for (k in 1:length(cell_number_list)){
cell_number=cell_number_list[[k]]  

for (j in 1:length(pbmc.list)){

pbmc=pbmc.list[[j]]
celltype=names(pbmc.list)[j]

Idents(pbmc) <- "projid"
pbmc=subset(pbmc, downsample = cell_number)
Idents(pbmc) <- "cell_type_high_resolution"

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("G:/PFC_429_final_Spring2021/Vanshika")
CRscores=read.csv("PFC427_CDR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

#muscat DE analysis
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb2))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

res <- pbDS(pb2,design = design,coef=as.list(2:4))

#Define output directory
output_directory=paste0("G:/PFC_429_final_Spring2021/Vanshika/muscat/Revision_Dec2022/Downsample_Exc/Results/Results_",cell_number,"_cells/gpath")
setwd(output_directory)

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_","gpath",".csv"))
}
}
}
##############################################

