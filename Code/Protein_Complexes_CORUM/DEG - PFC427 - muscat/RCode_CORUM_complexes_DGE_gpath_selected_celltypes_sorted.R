library(dplyr)
library(qpcR)

# coreComplexes

corecomplexes=read.csv("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Selected_coreComplexes/coreComplexes_Selected.csv")
corecomplexes=corecomplexes[corecomplexes$ComplexID %in% c("282","2839","165","166","3197","786","5464","7313","7567","6529","6166","7459","7220","6994"),]
corecomplexes=corecomplexes[,c("ComplexName", "subunits.Gene.name.")]

data_split=split(corecomplexes,corecomplexes$ComplexName)
data_split=lapply(data_split,"[", ,"subunits.Gene.name.")
data_split=lapply(data_split, strsplit,";", fixed = TRUE,"[[")
data_split=lapply(data_split, unlist,"[[")

Marker_genes_list_1=data_split
Marker_genes_list_1=Marker_genes_list_1[lengths(Marker_genes_list_1)>0]


# gpath DGE

setwd("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

list_celltypes=list.files()
#list_celltypes=list.files(,pattern = "Exc|Ast_|Mic_|Oli_|OPC")
list_celltypes=list_celltypes[c(4,7:11,14,18,52,53)]

for(i in 1:length(Marker_genes_list_1)){
  print(names(Marker_genes_list_1[i]))
  
  complex=as.data.frame(sort(Marker_genes_list_1[[i]]))
  colnames(complex)="gene"
  
  complex_2=data.frame()
  
  for (j in 1:length(list_celltypes)){
    celltype=read.csv(list_celltypes[j])
    celltype=left_join(complex,celltype,by="gene")
    
    celltype$logp=-log10(celltype$p_adj.loc)
    celltype$logp[celltype$logp<0]<-0
    celltype$signed_logp=celltype$logp*(celltype$logFC/(abs(celltype$logFC)))
    celltype=celltype[,c("gene","signed_logp")]
    
    celltype=as.data.frame(t(celltype))
    colnames(celltype)=celltype[1,]
    celltype=celltype[-1,drop=FALSE,]
    
    celltype_name=gsub("_gpath.csv","",list_celltypes[j])
    celltype=cbind(celltype_name,celltype)
      
    if(ncol(complex_2)>0) complex_2=rbind(complex_2,celltype)
    if(ncol(complex_2)==0) complex_2=celltype

  }
  
  
  write.csv(complex_2,file=paste0("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Selected_coreComplexes/muscat/",names(Marker_genes_list_1[i]),"_gpath_DEG_selected_cell_types.csv"))
  assign(names(Marker_genes_list_1[i]),complex_2)
  
}

