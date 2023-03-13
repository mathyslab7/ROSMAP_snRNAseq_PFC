library(data.table)
library(dplyr)
library(transite)
library(harmonicmeanp)
library(metap)

setwd("G:/PFC_429_final_Spring2021/Ghada/CORUM/coreComplexes/glm/gpath_BH/Exc")

adj_pval=read.csv("BH_pval_glm_PFC427_Exc_CORUM_gpath_BH.csv",row.names = 1)
reg_coeff=read.csv("Regression_coefficient_glm_PFC427_Exc_CORUM_gpath_BH.csv",row.names = 1)

adj_pval=adj_pval[,1:12]
reg_coeff=reg_coeff[,1:12]

adj_pval=na.omit(adj_pval)

#####################
# metap
# sumlog
aggregate_pval=apply(adj_pval,1,sumlog) 
#aggregate_pval=apply(adj_pval,1,sump) 
#aggregate_pval=apply(adj_pval,1,meanp)

agg_pval=data.frame()
for(i in 1:length(aggregate_pval)){
  pval=as.data.frame(aggregate_pval[[i]][["p"]])
  rownames(pval)=names(aggregate_pval[i])
  colnames(pval)="agg_pval"
  
  if(length(agg_pval)>0) agg_pval=rbind(agg_pval,pval)
  if(length(agg_pval)==0) agg_pval=pval
  
}
#####################

agg_pval=cbind(ComplexID=rownames(agg_pval),agg_pval)

reg_coeff$avg_reg_coeff=rowSums(reg_coeff)/ncol(reg_coeff)
avg_reg_coeff=reg_coeff[,"avg_reg_coeff",drop=FALSE]
avg_reg_coeff=cbind(ComplexID=rownames(avg_reg_coeff),avg_reg_coeff)

Exc=left_join(avg_reg_coeff,agg_pval,by="ComplexID")
Exc$logp=-log10(Exc$agg_pval)
Exc$assoc_score=Exc$logp * (Exc$avg_reg_coeff/(abs(Exc$avg_reg_coeff)))

Exc$ComplexID=substr(Exc$ComplexID,2,nchar(Exc$ComplexID))
Exc$ComplexID=substr(Exc$ComplexID,1,nchar(Exc$ComplexID)-1)


# Complex name 

complex_name=read.csv("G:/PFC_429_final_Spring2021/Ghada/CORUM/coreComplexes/coreComplexes.csv")
complex_name=complex_name[complex_name$Organism %like% "Human",]
colnames(complex_name)[1]="ComplexID"
complex_name=complex_name[,1:2]
complex_name$ComplexID=as.character(complex_name$ComplexID)

Exc=left_join(Exc,complex_name, by="ComplexID")

setwd("G:/PFC_429_final_Spring2021/Ghada/CORUM/coreComplexes/glm/gpath_BH/Exc")
write.csv(Exc,file="Exc_CORUM_gpath_glm_agg_pval_avg_reg_coeff_metap_sumlog_12.csv")


