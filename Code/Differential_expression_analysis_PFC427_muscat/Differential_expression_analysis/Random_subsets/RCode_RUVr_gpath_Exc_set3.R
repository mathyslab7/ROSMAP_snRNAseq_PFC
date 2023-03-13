library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)

##########################################################################################################
setwd("G:/PFC_429_final_Spring2021/Vanshika")
random_subsets=read.csv("PFC427_projids_random_subsets_v4.csv")
random_subsets = random_subsets %>% mutate(random_subset2 = case_when(grepl("random_subset_1",random_subset)~"random_subset_1",
                                                                      grepl("random_subset_2",random_subset)~"random_subset_1",
                                                                      grepl("random_subset_3",random_subset)~"random_subset_2",
                                                                      grepl("random_subset_4",random_subset)~"random_subset_2"))

random_subsets %>% mutate_if(is.numeric, as.character) -> random_subsets
subset_list=split(random_subsets,random_subsets$random_subset2)
##########################################################################################################


setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set3_final.rds")

#Split into individual cell types
Idents(pbmc)<-"cell_type_high_resolution"
pbmc.list <- SplitObject(pbmc, split.by = "cell_type_high_resolution")

for (j in 1:length(pbmc.list)){

pbmc=pbmc.list[[j]]
celltype=names(pbmc.list)[j]

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

#Loop through subsets##############
for (k in 1:length(subset_list)){
  group=names(subset_list)[k]
  formula = ~gpath+pmi+age_death+msex
  cd = as.data.frame(colData(pb))
  cd[,"projid"]=rownames(cd)
  cd=left_join(subset_list[[k]],cd,by="projid")
  rownames(cd)=cd[,"projid"]
  cd2=cd[,c("gpath","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  pb2=pb[,rownames(design)]
  
  ##################################
  library(RUVSeq)
  library(zebrafishRNASeq)
  
  counts=pb2@assays@data@listData[[celltype]]
  
  genes <- rownames(counts)
  gpath=pb2@colData@listData[["gpath"]]
  set <- newSeqExpressionSet(as.matrix(counts),
                             phenoData = data.frame(gpath, row.names=pb2@colData@rownames))
  set
  
  design <- model.matrix(~gpath, data=pData(set))
  
  y <- DGEList(counts=counts(set))
  y <- calcNormFactors(y, method="TMM")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  
  set4 <- RUVr(set, genes, k=10, res)
  #pData(set4)
  w=pData(set4)
  pb2$w_1 <- w$W_1
  pb2$w_2 <- w$W_2
  pb2$w_3 <- w$W_3
  pb2$w_4 <- w$W_4
  pb2$w_5 <- w$W_5
  pb2$w_6 <- w$W_6
  pb2$w_7 <- w$W_7
  pb2$w_8 <- w$W_8
  pb2$w_9 <- w$W_9
  pb2$w_10 <- w$W_10
  
  ##################################
  
  #muscat DE analysis
  formula = ~gpath+pmi+age_death+msex+w_1+w_2+w_3+w_4+w_5+w_6+w_7+w_8+w_9+w_10
  cd = as.data.frame(colData(pb2))
  cd2=cd[,c("gpath","pmi","age_death","msex","w_1","w_2","w_3","w_4","w_5","w_6","w_7","w_8","w_9","w_10")]
  design=model.matrix(formula,cd2)
  
  res <- pbDS(pb2,design = design,coef=as.list(2:4))
  
  #Define output directory
  working_directory=paste0("G:/PFC_429_final_Spring2021/Vanshika/muscat/Revision_Dec2022/RUVr_random_subsets/Two_subsets/Results/gpath/",group)
  setwd(working_directory)
  
  #Indicate name of the variable analyzed
  variable="gpath"
  
  tbl <- res$table[[1]]
  names(tbl)
  
  for (i in 1:length(tbl)){
    name=names(tbl)[i]
    name=gsub("/","_",name)
    k1 <- tbl[[i]]
    write.csv(k1,file=paste0(name,"_","gpath",".csv"))
  }
}
}
##############################################