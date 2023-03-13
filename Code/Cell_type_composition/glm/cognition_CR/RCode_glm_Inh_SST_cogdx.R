library(dplyr)

variables=c(
  "projid",
  "cogdx",
  "AD",
  "msex",                        
  "age_death",
  "pmi")

cell_types=c("projid","Inh.SST")
setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12")

data_input=read.csv("PFC427_meta_number_of_cells_ordered_proportion_analyzed_metadata.csv")
data_input[,"Fib"]=data_input[,"Fib.FLRT2"]+data_input[,"Fib.SLC4A4"]

celltypes=data_input[,cell_types]
variables_info=data_input[,variables]

data=left_join(celltypes,variables_info,by="projid")
data=subset(data,AD=="yes")
data=subset(data,cogdx %in% c(1,4))

#Add binary coding for cogdx
data = data %>% mutate(cogdx_binary =
                                 case_when(cogdx == "4" ~ 0, 
                                           cogdx == "1" ~ 1)
)

#########################################################################
#Cell type
var="Inh.SST"
#Variable
var2="cogdx_binary"


#formula
fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death","+","msex"))
# fit glm model
fit <- glm(fmla, data=data,family = quasibinomial())
confint_table=confint(fit)
# get coefficents of fit
cfit <- coef(summary(fit))
    
# create data frame
df <- data.frame(cell_type = var, variable=var2, Estimate = cfit[var2,'Estimate'],CI_2.5=confint_table[var2,'2.5 %'],CI_97.5=confint_table[var2,'97.5 %'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'],
                     Deviance = deviance(fit), stringsAsFactors = F)
    
results <- df

setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Revision_Dec2022/Cell_fractions_glm")
write.csv(results,file="Results_glm_Inh_SST.csv")

