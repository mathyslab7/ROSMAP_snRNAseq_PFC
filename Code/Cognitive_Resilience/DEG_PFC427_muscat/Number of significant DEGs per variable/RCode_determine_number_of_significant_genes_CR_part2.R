library(Seurat)
library(muscat)
library(dplyr)
library(data.table)

setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results")

#Select variables
variable_list=list(
  "tangles_CR_score",
  "gpath_CDR_score",
  "nft_CDR_score",
  "plaq_n_CDR_score",
  "tangles_CDR_score"
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
print(variable)

pattern=paste0("*",variable,".RData")
list.filenames=list.files()
list.filenames=list.filenames[list.filenames %like% pattern]

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
  print(variable)
  
  pattern=paste0("*",variable,".RData")
  list.filenames=list.files()
  list.filenames=list.filenames[list.filenames %like% pattern]
  
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


setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Number_of_significant_DEGs")
write.csv(number_of_DEGs,file = "DGE_Comprehensive_number_of_significant_genes_CR_part2.csv")

######################################

setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Number_of_significant_DEGs")
number_of_DEGs_1=read.csv("DGE_Comprehensive_number_of_significant_genes_CR_part1.csv",row.names = 1)
number_of_DEGs_2=number_of_DEGs

number_of_DEGs=left_join(number_of_DEGs_1,number_of_DEGs_2,by="cell_type_high_resolution")
write.csv(number_of_DEGs,file = "DGE_Comprehensive_number_of_significant_genes_CR.csv")

order=read.csv("G:/PFC_429_final_Spring2021/Ghada/Celltypes Annotations and Numbers/PFC_celltypes_order.csv")
colnames(order)[1]="cell_type_high_resolution_orig"
order$cell_type_high_resolution=make.names(order$cell_type_high_resolution_orig)

number_of_DEGs=left_join(order,number_of_DEGs,by="cell_type_high_resolution")
rownames(number_of_DEGs)=number_of_DEGs[,1]
number_of_DEGs=number_of_DEGs[,-c(1,2)]

write.csv(number_of_DEGs,file = "DGE_Comprehensive_number_of_significant_genes_CR_ordered.csv")


