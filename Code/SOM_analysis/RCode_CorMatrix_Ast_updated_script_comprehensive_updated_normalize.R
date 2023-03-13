#Run updated code

setwd("F:/PFC_429_final_Spring2021/Data")
library(Seurat)
library(dplyr)
pbmc=readRDS("Ast_integrated_batch_5000_leiden_predicted_final.rds")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("F:/PFC_429_final_Spring2021/Vanshika")

metadata=read.csv("Combined_metadata.csv")
metadata=metadata[-1]
metadata=metadata[-2]
metadata=metadata[-2]
#Make sure 39 is the right column
cols=colnames(metadata)
cols
metadata=metadata[-39]
metadatanew=metadata[,c(1,11,13:15,16:27,29,31,54:57,61:62,67:69,72:76,80,83:86,90,93,108,112,114,115,117:125)]
metadata=metadatanew

################################################################################
#Add additional CR scores and remove age_death CR score
setwd("F:/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
CRscores=CRscores[,c("projid","plaq_n_CR_score","nft_CR_score")]

library(dplyr)
New_metadata=left_join(metadata,CRscores,by='projid')
New_metadata=New_metadata[-53]
metadata=New_metadata

################################################################################
#Adjust variables
metadata=metadata %>% mutate(apoe_genotype = case_when(apoe_genotype == '22' ~ 0,
                                                       apoe_genotype == '23' ~ 0,
                                                       apoe_genotype == '33' ~ 0,
                                                       apoe_genotype == '24' ~ 1,
                                                       apoe_genotype == '34' ~ 1,
                                                       apoe_genotype == '44' ~ 1))

metadata=metadata %>% mutate(cogdx = na_if(cogdx, "3"))
metadata=metadata %>% mutate(cogdx = na_if(cogdx, "5"))
metadata=metadata %>% mutate(cogdx = na_if(cogdx, "6"))

metadata=metadata %>% mutate(cogdx_stroke = na_if(cogdx_stroke, "9"))

metadata=metadata %>% mutate(dxpark = na_if(dxpark, "9"))

metadata=metadata %>% mutate(dlbdx = na_if(dlbdx, "1"))
metadata=metadata %>% mutate(dlbdx = na_if(dlbdx, "2"))
################################################################################

object.list=SplitObject(pbmc,split.by = "cell_type_high_resolution")
Average_expression_list<-list()
for (i in 1:length(object.list)){
  name=names(object.list)[i]
  average=AverageExpression(object.list[[i]],slot = "data",assays = "RNA",group.by="projid")
  RNA_Average = average$RNA
  RNA_Average_Dataframe = as.data.frame(RNA_Average)
  Average_expression_list[[i]]=RNA_Average_Dataframe
  names(Average_expression_list)[i]<-name
}

Df_list<-list()
for (i in 1:length(Average_expression_list)){
  name=names(Average_expression_list)[i]
  average=Average_expression_list[[i]]
  average_t=t(average)
  expression_df=as.data.frame(average_t)
  df2 <- cbind(projid = rownames(average_t), average_t)
  df3=as.data.frame(df2)
  df4 <- df3[c("projid")]
  df4$projid <- as.integer(as.character(df4$projid))
  df5=left_join(df4,metadata,by="projid")
  df6=df5[-1]
  
  CorMatrix=cor(df6,expression_df,method = 'spearman',use="pairwise.complete.obs")
  df10=as.data.frame(CorMatrix)
  Df_list[[i]]=df10
  names(Df_list)[i]=name
}

Df_list2<-list()
for (i in 1:length(Df_list)){
  name=names(Df_list)[i]
  df11=as.data.frame(Df_list[[i]])
  df11_t=t(df11)
  colnames(df11_t) <- paste(name, colnames(df11_t), sep = "_")
  df12=as.data.frame(df11_t)
  Df_list2[[i]]=df12
  names(Df_list2)[i]=name
}

test=as.data.frame(Df_list2)
setwd("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021")
write.csv(test,file="CorMatrix_Ast_updated_script_comprehensive.csv")