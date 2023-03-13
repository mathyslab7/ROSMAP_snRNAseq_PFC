library(dplyr)
setwd("D:/PFC_427_WS/Vanshika")

celltypes=read.csv("celltype_order_Mic_Ast_namesR.csv",row.names = 1)

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/MeanModuleScores/GLM_results")
data=read.csv("MeanModuleScores_complex165_glm_results.csv",row.names = 1)
rownames(data)=data$variable


datat=t(data)
datat=as.data.frame(datat)
datat=datat[-1,]

datat$celltypeR=rownames(datat)

data_new=left_join(celltypes,datat,by="celltypeR")

write.csv(data_new,file="MeanModuleScores_complex165_glm_results_ordered.csv")