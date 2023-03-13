setwd("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021")
exset3=read.csv(file="CorMatrix_Exc_raw_set3_final_updated_script_comprehensive.csv")
exset2=read.csv(file="CorMatrix_Exc_raw_set2_final_updated_script_comprehensive.csv")
exset1=read.csv(file="CorMatrix_Exc_raw_set1_final_updated_script_comprehensive.csv")
opc=read.csv(file="CorMatrix_opc_updated_script_comprehensive.csv")
oli=read.csv(file="CorMatrix_Oli_updated_script_comprehensive.csv")
immune=read.csv(file="CorMatrix_immune_updated_script_comprehensive.csv")
inh=read.csv(file="CorMatrix_inh_updated_script_comprehensive.csv")
ast=read.csv(file="CorMatrix_Ast_updated_script_comprehensive.csv")
vasc=read.csv(file="CorMatrix_Vasc_updated_script_comprehensive.csv")
library(dplyr)
combine=left_join(exset1,exset2,by="X")
combine_2=left_join(combine,exset3,by="X")
combine_3=left_join(combine_2,opc,by="X")
combine_4=left_join(combine_3,oli,by="X")
combine_5=left_join(combine_4,immune,by="X")
combine_6=left_join(combine_5,inh,by="X")
combine_7=left_join(combine_6,ast,by="X")
combine_all=left_join(combine_7,vasc,by="X")
rownames(combine_all)=combine_all$X
combine_all=combine_all[-1]
#output
#write.csv(combine_all,file="combine_all_dplyr.csv")

CorMatrix_new=combine_all[-which(rowMeans(is.na(combine_all)) > 0.5), ]

write.csv(CorMatrix_new,file="CorMatrix_dplyr_comprehensive.csv")