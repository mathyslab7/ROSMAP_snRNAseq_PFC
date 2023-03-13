library(Seurat)
library(muscat)
library(dplyr)

setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results")

#Select variables
variable_list=list(
  "cogdx_AD_yes",
  "gpath_CR_score",
  "nft_CR_score",
  "plaq_n_CR_score"
)

for (i in 1:length(variable_list)){
    if(grepl("cogdx_AD_yes",variable_list[[i]],fixed=TRUE)){
      variable="cogdx"
    } else {
      variable=variable_list[[i]]
    }  
    print(variable)
}

###########################################
variable=variable_list[[1]]

pattern=paste0("*",variable,".RData")
list.filenames=list.files(pattern=pattern)

#number_of_DEGs=data.frame()

load(list.filenames[1])
# access results table for 1st comparison
tbl <- res$table[[1]]
for (i in 1:length(tbl)) {
  tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
  tbl[[i]]=nrow(tbl[[i]])
}
number_of_DEGs<-as.data.frame(tbl)

for (i in 2:length(list.filenames)){
  load(list.filenames[i])
  # access results table for 1st comparison
  tbl <- res$table[[1]]
  for (i in 1:length(tbl)) {
    tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
    tbl[[i]]=nrow(tbl[[i]])
  }
  DF<-as.data.frame(tbl)
  number_of_DEGs=cbind(number_of_DEGs,DF)
}
rownames(number_of_DEGs)<-variable_list[[1]]
number_of_DEGs=t(number_of_DEGs)
number_of_DEGs=cbind(cell_type_high_resolution=rownames(number_of_DEGs),number_of_DEGs)
number_of_DEGs=as.data.frame(number_of_DEGs)


for (j in 2:length(variable_list)) {
  
  #Select variable (here cogdx Workspaces)
  variable=variable_list[[j]]
  
  pattern=paste0("*",variable,".RData")
  list.filenames=list.files(pattern=pattern)
  
  #number_of_DEGs=data.frame()
  
  load(list.filenames[1])
  # access results table for 1st comparison
  tbl <- res$table[[1]]
  for (i in 1:length(tbl)) {
    tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
    tbl[[i]]=nrow(tbl[[i]])
  }
  number_of_DEGs2<-as.data.frame(tbl)
  
  for (i in 2:length(list.filenames)){
    load(list.filenames[i])
    # access results table for 1st comparison
    tbl <- res$table[[1]]
    for (i in 1:length(tbl)) {
      tbl[[i]]=subset(tbl[[i]],p_adj.loc < 0.05)
      tbl[[i]]=nrow(tbl[[i]])
    }
    DF<-as.data.frame(tbl)
    number_of_DEGs2=cbind(number_of_DEGs2,DF)
  }
  rownames(number_of_DEGs2)<-variable_list[[j]]
  number_of_DEGs2=t(number_of_DEGs2)
  number_of_DEGs2=cbind(cell_type_high_resolution=rownames(number_of_DEGs2),number_of_DEGs2)
  number_of_DEGs2=as.data.frame(number_of_DEGs2)
  
  number_of_DEGs=full_join(number_of_DEGs,number_of_DEGs2,by="cell_type_high_resolution")
}


setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Number_of_significant_DEGs")
write.csv(number_of_DEGs,file = "DGE_Comprehensive_number_of_significant_genes_CR_part1.csv")
