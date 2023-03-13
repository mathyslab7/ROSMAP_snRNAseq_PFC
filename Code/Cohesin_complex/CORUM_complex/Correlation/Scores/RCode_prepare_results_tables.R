library(dplyr)
setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores")

filelist=list.files(pattern="*score.csv")

############################################################################
for (i in 1:length(filelist)){
  name=tools::file_path_sans_ext(filelist[[i]])
  name=gsub("_score","",name)
  data=read.csv(filelist[[i]])
  colnames(data)[colnames(data) == "X"] <- "gene"
  data2=data[,c("gene","X1651")]
  colnames(data2)[colnames(data2) == "X1651"] <- name
  if (i==1){
    results_table=data2
  }
  if (i>1){
    results_table=full_join(results_table,data2,by="gene")
  }
}

write.csv(results_table,file="Results_table_scores_complex165.csv")
############################################################################

############################################################################
for (i in 1:length(filelist)){
  name=tools::file_path_sans_ext(filelist[[i]])
  name=gsub("_score","",name)
  data=read.csv(filelist[[i]])
  colnames(data)[colnames(data) == "X"] <- "gene"
  data2=data[,c("gene","X1661")]
  colnames(data2)[colnames(data2) == "X1661"] <- name
  if (i==1){
    results_table=data2
  }
  if (i>1){
    results_table=full_join(results_table,data2,by="gene")
  }
}

write.csv(results_table,file="Results_table_scores_complex166.csv")
############################################################################

############################################################################
for (i in 1:length(filelist)){
  name=tools::file_path_sans_ext(filelist[[i]])
  name=gsub("_score","",name)
  data=read.csv(filelist[[i]])
  colnames(data)[colnames(data) == "X"] <- "gene"
  data2=data[,c("gene","X54641")]
  colnames(data2)[colnames(data2) == "X54641"] <- name
  if (i==1){
    results_table=data2
  }
  if (i>1){
    results_table=full_join(results_table,data2,by="gene")
  }
}

write.csv(results_table,file="Results_table_scores_complex5464.csv")
############################################################################