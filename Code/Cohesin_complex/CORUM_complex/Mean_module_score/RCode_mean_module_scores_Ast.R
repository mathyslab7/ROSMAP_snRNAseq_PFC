library(dplyr)
setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Cell_numbers_fractions")
MainFrame=read.csv("PFC427_projids.csv")


setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Module_scores")

data=read.csv("Ast_integrated_batch_5000_leiden_predicted_final_metadata.csv",row.names=1)
colnames(data)

data=data[,c("projid","cell_type_high_resolution","X1651")]

data_split=split(data, data$cell_type_high_resolution)

for (i in 1:length(data_split)) {
  name=names(data_split)[i]
  data2=data_split[[i]]
  data_split2=split(data2, data2$projid)
  for (j in 1:length(data_split2)) {
    name2=names(data_split2)[j]
    data3=data_split2[[j]]
    data4=data3[,"X1651",drop=FALSE]
    data5=colMeans(data4)
    TempDf=data.frame(data5)
    colnames(TempDf)<-name
    rownames(TempDf)<-name2
    name2=as.numeric(name2)
    TempDf[,"projid"]<-name2
    if (j==1){
      TempResults=TempDf
    }
    if (j>1){
      TempResults=rbind(TempResults,TempDf)
    }
  }
  MainFrame=left_join(MainFrame,TempResults,by="projid")
}

#Major cell type (Ast)

data_split3=split(data, data$projid)

for (k in 1:length(data_split3)) {
    name2=names(data_split3)[k]
    data3=data_split3[[k]]
    data4=data3[,"X1651",drop=FALSE]
    data5=colMeans(data4)
    TempDf=data.frame(data5)
    colnames(TempDf)<-"Ast"
    rownames(TempDf)<-name2
    name2=as.numeric(name2)
    TempDf[,"projid"]<-name2
    if (k==1){
      TempResults=TempDf
    }
    if (k>1){
      TempResults=rbind(TempResults,TempDf)
    }
  }
MainFrame=left_join(MainFrame,TempResults,by="projid")

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/MeanModuleScores")
write.csv(MainFrame,file="MeanModuleScores_complex165_Ast.csv")


