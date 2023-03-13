library(dplyr)

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/MeanModuleScores")
Vasc=read.csv(file="MeanModuleScores_complex165_Vasc.csv")
Vasc=subset (Vasc, select = -X)

OPC=read.csv(file="MeanModuleScores_complex165_OPC.csv")
OPC=subset (OPC, select = -X)

Oli=read.csv(file="MeanModuleScores_complex165_Oli.csv")
Oli=subset (Oli, select = -X)

Mic=read.csv(file="MeanModuleScores_complex165_Mic.csv")
Mic=subset (Mic, select = -X)

Inh=read.csv(file="MeanModuleScores_complex165_Inh.csv")
Inh=subset (Inh, select = -X)

Exc3=read.csv(file="MeanModuleScores_complex165_Exc3.csv")
Exc3=subset (Exc3, select = -X)

Exc2=read.csv(file="MeanModuleScores_complex165_Exc2.csv")
Exc2=subset (Exc2, select = -X)

Exc1=read.csv(file="MeanModuleScores_complex165_Exc1.csv")
Exc1=subset (Exc1, select = -X)

Ast=read.csv(file="MeanModuleScores_complex165_Ast.csv")
Ast=subset (Ast, select = -X)

combined=left_join(Ast,Exc1,by="projid")
combined2=left_join(combined,Exc2,by="projid")
combined3=left_join(combined2,Exc3,by="projid")
combined4=left_join(combined3,Inh,by="projid")
combined5=left_join(combined4,Mic,by="projid")
combined6=left_join(combined5,Oli,by="projid")
combined7=left_join(combined6,OPC,by="projid")
combined8=left_join(combined7,Vasc,by="projid")

write.csv(combined8,file="MeanModuleScores_complex165.csv")

