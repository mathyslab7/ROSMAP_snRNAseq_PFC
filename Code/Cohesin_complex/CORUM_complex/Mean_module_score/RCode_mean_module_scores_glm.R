library(dplyr)

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/MeanModuleScores")
data=read.csv(file="MeanModuleScores_complex165.csv")
data=subset (data, select = -X)

setwd("D:/Projects/Large_scale_snRNA_seq/PFC_429_final/Metadata/Confirmation_Cell_numbers")
metadata=read.csv(file="Combined_metadata.csv")

data_table=left_join(data,metadata,by="projid")

#Replacing values that do not make sense to be incorporated in the analysis
data_table <- data_table %>%                               # Replacing values
  mutate(dxpark = replace(dxpark, dxpark == 9, NA))

data_table <- data_table %>%                               # Replacing values
  mutate(cogdx_stroke = replace(cogdx_stroke, cogdx_stroke == 9, NA))

data_table <- data_table %>%                               # Replacing values
  mutate(dlbdx = replace(dlbdx, dlbdx == 1, NA))

data_table <- data_table %>%                               # Replacing values
  mutate(dlbdx = replace(dlbdx, dlbdx == 2, NA))

data_table <- data_table %>%                               # Replacing values
  mutate(Apoe_e4 = replace(Apoe_e4, Apoe_e4 == "unknown", NA))

data_table <- data_table %>%                               # Replacing values
  mutate(Apoe_e4 = replace(Apoe_e4, Apoe_e4 == "yes", 1))

data_table <- data_table %>%                               # Replacing values
  mutate(Apoe_e4 = replace(Apoe_e4, Apoe_e4 == "no", 0))

#Adjust categorical variables such that disease is encoded by 1 and control is encoded by 0
data_table <- data_table %>%                               # Replacing values
  mutate(cogdx_stroke = replace(cogdx_stroke, cogdx_stroke == 1, 1))

data_table <- data_table %>%                               # Replacing values
  mutate(cogdx_stroke = replace(cogdx_stroke, cogdx_stroke == 2, 0))

data_table <- data_table %>%                               # Replacing values
  mutate(dxpark = replace(dxpark, dxpark == 1, 1))

data_table <- data_table %>%                               # Replacing values
  mutate(dxpark = replace(dxpark, dxpark == 2, 0))

data_table=subset(data_table,select=c(1:57,73:74,90,114:117,121,122,127:129,134,144:146,153,168,174:175,177:183,195,88,143))

# create data frame to store results
results <- data.frame()

# loop through the scales and each variable
for(var in names(data_table)[c(2:57)]){
  for(var2 in names(data_table)[c(58:85)]){
    # dynamically generate formula
    fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death"))
    # fit glm model
    fit <- glm(fmla, data=data_table)
    
    # get coefficents of fit
    cfit <- coef(summary(fit))
    
   if (var2=="Apoe_e4"){
      var2="Apoe_e41"
    }
    
    # create temporary data frame
    df <- data.frame(cell_type = var, variable=var2, Estimate = cfit[var2,'Estimate'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'],
                     Deviance = deviance(fit), stringsAsFactors = F)
    
    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
  }
}

results[,"p_adj"]=p.adjust(results[,"p.value"],method ="BH")
results[,"log10_p_adj"]=-log10(results[,"p_adj"])
results[,"score"]=results[,"log10_p_adj"]*(results[,"Estimate"]/abs(results[,"Estimate"]))

data_split=split(results,results$cell_type)
for (i in 1:length(data_split)){
  name=names(data_split)[i]
  temp_data=data_split[[i]]
  temp_data=temp_data[,c("variable","score")]
  names(temp_data)[names(temp_data) == "score"] <- name
  if (i==1){
    Results_table=temp_data
  }
  if (i>1){
    Results_table=full_join(Results_table,temp_data,by="variable")
  }
}

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/MeanModuleScores/GLM_results")
write.csv(Results_table,file="MeanModuleScores_complex165_glm_results.csv")