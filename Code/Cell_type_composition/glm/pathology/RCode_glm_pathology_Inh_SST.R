library(dplyr)

variables=c(
  "projid",
  "gpath",
  "gpath_3neocort",                  
  "plaq_d_mf",
  "amyloid",
  "plaq_d",
  "plaq_n",
  "plaq_n_mf",
  "nft",
  "nft_mf",
  "tangles",
  "msex",                        
  "age_death",
  "pmi")

cell_types=c("projid",
             "Inh.CUX2.MSR1",                   
             "Inh.ENOX2.SPHKAP",
             "Inh.L3.5.SST.MAFB")
setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12")

data_input=read.csv("PFC427_meta_number_of_cells_ordered_proportion_analyzed_metadata.csv")
data_input[,"Fib"]=data_input[,"Fib.FLRT2"]+data_input[,"Fib.SLC4A4"]

celltypes=data_input[,cell_types]
variables_info=data_input[,variables]

data=left_join(celltypes,variables_info,by="projid")

#########################################################################
variable_list=c(
  "gpath",
  "gpath_3neocort",                  
  "plaq_d_mf",
  "amyloid",
  "plaq_d",
  "plaq_n",
  "plaq_n_mf",
  "nft",
  "nft_mf",
  "tangles")

cell_types_list=c("Inh.CUX2.MSR1",                   
                  "Inh.ENOX2.SPHKAP",
                  "Inh.L3.5.SST.MAFB")

for (i in 1:length(variable_list)){
  # create data frame to store results
  results <- data.frame()
for(var in cell_types_list){
  var2 = variable_list[i]
    # dynamically generate formula
    fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death","+","msex"))
    # fit glm model
    fit <- glm(fmla, data=data,family = quasibinomial())
    confint_table=confint(fit)
    # get coefficents of fit
    cfit <- coef(summary(fit))
    
    # create temporary data frame
    df <- data.frame(cell_type = var, variable=var2, Estimate = cfit[var2,'Estimate'],CI_2.5=confint_table[var2,'2.5 %'],CI_97.5=confint_table[var2,'97.5 %'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'],
                     Deviance = deviance(fit), stringsAsFactors = F)
    
    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
  }
  results[,"pval_adj"]=p.adjust(results[,"p.value"],method="BH")
  setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12/glm/Inh_SST_pathology")
  file_name=paste0("Inh_subclass_",var2,".csv")
  write.csv(results,file=file_name)
}
