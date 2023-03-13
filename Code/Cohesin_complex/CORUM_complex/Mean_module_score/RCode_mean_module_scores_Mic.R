library(dplyr)
setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Cell_numbers_fractions")
MainFrame=read.csv("PFC427_projids.csv")


setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Module_scores")

data=read.csv("Immune_integrated_batch_5000_predicted_final_metadata.csv",row.names=1)
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

#Major cell type (Mic)

data_2 = data %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh ",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast ",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli ",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC ",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic ",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))

data_2=subset(data_2,major_cell_type == "Mic")


data_split3=split(data_2, data_2$projid)

for (k in 1:length(data_split3)) {
    name2=names(data_split3)[k]
    data3=data_split3[[k]]
    data4=data3[,"X1651",drop=FALSE]
    data5=colMeans(data4)
    TempDf=data.frame(data5)
    colnames(TempDf)<-"Mic"
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
write.csv(MainFrame,file="MeanModuleScores_complex165_Mic.csv")


