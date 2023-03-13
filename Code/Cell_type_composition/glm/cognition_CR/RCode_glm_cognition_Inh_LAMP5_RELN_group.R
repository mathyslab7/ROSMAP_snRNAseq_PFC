library(dplyr)

variables=c(
  "projid",
  "cogn_global_lv",
  "cogn_ep_lv",                      
  "cogn_po_lv",
  "cogn_ps_lv",
  "cogn_se_lv",
  "cogn_wo_lv",   
  "cogng_random_slope",
  "cognep_random_slope",
  "cognpo_random_slope",
  "cognps_random_slope",             
  "cognse_random_slope",
  "cognwo_random_slope",
  "msex",                        
  "age_death",
  "pmi")

cell_types=c("projid",
             "Inh.LAMP5.RELN",
             "Inh.SORCS1.TTN",                                    
             "Inh.PTPRK.FAM19A1")
setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12")

data_input=read.csv("PFC427_meta_number_of_cells_ordered_proportion_analyzed_metadata.csv")
data_input[,"Fib"]=data_input[,"Fib.FLRT2"]+data_input[,"Fib.SLC4A4"]

celltypes=data_input[,cell_types]
variables_info=data_input[,variables]

data=left_join(celltypes,variables_info,by="projid")

#########################################################################
variable_list=c(
  "cogn_global_lv",
  "cogn_ep_lv",                      
  "cogn_po_lv",
  "cogn_ps_lv",
  "cogn_se_lv",
  "cogn_wo_lv",   
  "cogng_random_slope",
  "cognep_random_slope",
  "cognpo_random_slope",
  "cognps_random_slope",             
  "cognse_random_slope",
  "cognwo_random_slope")

cell_types_list=c("Inh.LAMP5.RELN",
                  "Inh.SORCS1.TTN",                                    
                  "Inh.PTPRK.FAM19A1")

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
  setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Metadata_from_Seurat_objects/PFC427_Exc12/glm/Inh_LAMP5_RELN_group_cognition")
  file_name=paste0("Inh_subclass_",var2,".csv")
  write.csv(results,file=file_name)
}
