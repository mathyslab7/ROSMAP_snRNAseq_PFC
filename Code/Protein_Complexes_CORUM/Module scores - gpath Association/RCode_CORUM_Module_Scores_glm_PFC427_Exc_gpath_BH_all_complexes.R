library(dplyr)
library(data.table)

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
data=readRDS("Exc_raw_set1_final_metadata_relevant.rds")

cell_type_split=split(data,data$cell_type_high_resolution)


table_pval_Exc1=data.frame()
table_coeff_Exc1=data.frame()

for(cell in 1:length(cell_type_split)){
  
  data_split=split(cell_type_split[[cell]], cell_type_split[[cell]]$projid)
  
  #First element
  name=names(data_split)[1]
  data3=data_split[[1]]
  data4=data3[,c(3:2341)]
  data5=colMeans(data4)
  data6=as.data.frame(data5)
  colnames(data6)=name
  main_table=data6
  
  for (i in 2:length(data_split)) {
    name=names(data_split)[i]
    data3=data_split[[i]]
    data4=data3[,c(3:2341)]
    data5=colMeans(data4)
    data6=as.data.frame(data5)
    colnames(data6)=name
    main_table=cbind(main_table,data6)
  }
  
  main_table_new=t(main_table)
  main_table_new=as.data.frame(main_table_new)
  main_table_new2=cbind(projid=rownames(main_table_new),main_table_new)
  
  setwd("F:/PFC_429_final_Spring2021/Vanshika")
  metadata=read.csv("dataset_652_basic_04-23-2020.csv")
  metadata1 <- metadata %>%  
    mutate(projid = as.character(projid))
  
  main_table_anno=left_join(main_table_new2,metadata1,by="projid")
  
  data_table=select(main_table_anno,1:2340,"age_death","pmi","gpath")
  

  
  
  #create data frame to store results
  results <- data.frame()
  
  results_coeff_final <- data.frame()
  results_pval_final <- data.frame()
  
  #loop through module scores and variables
  for(var in names(data_table)[c(2:2340)]){
    results_coeff <- data.frame()
    results_pval <- data.frame()
    for(var2 in names(data_table)[2343]){
      # dynamically generate formula
      fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death"))
      # fit glm model
      fit <- glm(fmla, data=data_table)
      
      # get coefficents of fit
      cfit <- coef(summary(fit))
      
      # create temporary data frame
      df <- data.frame(module_score = var, variable=var2, coefficient = cfit[var2,'Estimate'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
      df_coeff <- data.frame(coefficient = cfit[var2,'Estimate'], stringsAsFactors = F)
      rownames(df_coeff)=var
      colnames(df_coeff)=var2
      df_pval <- data.frame(p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
      rownames(df_pval)=var
      colnames(df_pval)=var2
      
      # bind rows of temporary data frame to the results data frame
      results <- rbind(results, df)
      
      #################################################
      if (nrow(results_coeff)>0) {
        results_coeff <- cbind(results_coeff, df_coeff)
      }
      if (nrow(results_coeff)==0) {
        results_coeff <- df_coeff
      }
      
      cell_type_name=names(cell_type_split[cell])
      colnames(results_coeff)=cell_type_name
      #################################################
      
      #################################################
      if (nrow(results_pval)>0) {
        results_pval <- cbind(results_pval, df_pval)
      }
      if (nrow(results_pval)==0) {
        results_pval <- df_pval
      }
      
      cell_type_name=names(cell_type_split[cell])
      colnames(results_pval)=cell_type_name
      #################################################
    }
    #################################################
    if (ncol(results_coeff_final)>0) {
      results_coeff_final <- rbind(results_coeff_final, results_coeff)
    }
    if (ncol(results_coeff_final)==0) {
      results_coeff_final <- results_coeff
    }
    
    #################################################
    
    #################################################
    if (ncol(results_pval_final)>0) {
      results_pval_final <- rbind(results_pval_final, results_pval)
    }
    if (ncol(results_pval_final)==0) {
      results_pval_final <- results_pval
    }
    
    #################################################
  }
  
  # table_coeff_Exc1 (all cell types)
  if(length(table_coeff_Exc1)>0){
    table_coeff_Exc1=cbind(table_coeff_Exc1,results_coeff_final)
  }
  
  if(length(table_coeff_Exc1)==0){
    table_coeff_Exc1=results_coeff_final
  }
  
  
  # table_pval_Exc1 (all cell types)
  if(length(table_pval_Exc1)>0){
    table_pval_Exc1=cbind(table_pval_Exc1,results_pval_final)
  }
  
  if(length(table_pval_Exc1)==0){
    table_pval_Exc1=results_pval_final
  }
}



#################################################################################
#################################################################################

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
data=readRDS("Exc_raw_set2_final_metadata_relevant.rds")

cell_type_split=split(data,data$cell_type_high_resolution)


table_pval_Exc2=data.frame()
table_coeff_Exc2=data.frame()

for(cell in 1:length(cell_type_split)){

data_split=split(cell_type_split[[cell]], cell_type_split[[cell]]$projid)

#First element
name=names(data_split)[1]
data3=data_split[[1]]
data4=data3[,c(3:2341)]
data5=colMeans(data4)
data6=as.data.frame(data5)
colnames(data6)=name
main_table=data6

for (i in 2:length(data_split)) {
  name=names(data_split)[i]
  data3=data_split[[i]]
  data4=data3[,c(3:2341)]
  data5=colMeans(data4)
  data6=as.data.frame(data5)
  colnames(data6)=name
  main_table=cbind(main_table,data6)
}

main_table_new=t(main_table)
main_table_new=as.data.frame(main_table_new)
main_table_new2=cbind(projid=rownames(main_table_new),main_table_new)

setwd("F:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")
metadata1 <- metadata %>%  
  mutate(projid = as.character(projid))

main_table_anno=left_join(main_table_new2,metadata1,by="projid")

data_table=select(main_table_anno,1:2340,"age_death","pmi","gpath")




#create data frame to store results
results <- data.frame()

results_coeff_final <- data.frame()
results_pval_final <- data.frame()

#loop through module scores and variables
for(var in names(data_table)[c(2:2340)]){
  results_coeff <- data.frame()
  results_pval <- data.frame()
  for(var2 in names(data_table)[2343]){
    # dynamically generate formula
    fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death"))
    # fit glm model
    fit <- glm(fmla, data=data_table)
    
    # get coefficents of fit
    cfit <- coef(summary(fit))
    
    # create temporary data frame
    df <- data.frame(module_score = var, variable=var2, coefficient = cfit[var2,'Estimate'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
    df_coeff <- data.frame(coefficient = cfit[var2,'Estimate'], stringsAsFactors = F)
    rownames(df_coeff)=var
    colnames(df_coeff)=var2
    df_pval <- data.frame(p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
    rownames(df_pval)=var
    colnames(df_pval)=var2
    
    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
    
    #################################################
    if (nrow(results_coeff)>0) {
      results_coeff <- cbind(results_coeff, df_coeff)
    }
    if (nrow(results_coeff)==0) {
      results_coeff <- df_coeff
    }
    
    cell_type_name=names(cell_type_split[cell])
    colnames(results_coeff)=cell_type_name
    #################################################
    
    #################################################
    if (nrow(results_pval)>0) {
      results_pval <- cbind(results_pval, df_pval)
    }
    if (nrow(results_pval)==0) {
      results_pval <- df_pval
    }
    
    cell_type_name=names(cell_type_split[cell])
    colnames(results_pval)=cell_type_name
    #################################################
  }
  #################################################
  if (ncol(results_coeff_final)>0) {
    results_coeff_final <- rbind(results_coeff_final, results_coeff)
  }
  if (ncol(results_coeff_final)==0) {
    results_coeff_final <- results_coeff
  }
  
  #################################################
  
  #################################################
  if (ncol(results_pval_final)>0) {
    results_pval_final <- rbind(results_pval_final, results_pval)
  }
  if (ncol(results_pval_final)==0) {
    results_pval_final <- results_pval
  }
  
  #################################################
}

# table_coeff_Exc2 (all cell types)
if(length(table_coeff_Exc2)>0){
  table_coeff_Exc2=cbind(table_coeff_Exc2,results_coeff_final)
}

if(length(table_coeff_Exc2)==0){
  table_coeff_Exc2=results_coeff_final
}


# table_pval_Exc2 (all cell types)
if(length(table_pval_Exc2)>0){
  table_pval_Exc2=cbind(table_pval_Exc2,results_pval_final)
}

if(length(table_pval_Exc2)==0){
  table_pval_Exc2=results_pval_final
}
}


#################################################################################
#################################################################################

setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/Module_scores")
data=readRDS("Exc_raw_set3_final_metadata_relevant.rds")

cell_type_split=split(data,data$cell_type_high_resolution)


table_pval_Exc3=data.frame()
table_coeff_Exc3=data.frame()

for(cell in 1:length(cell_type_split)){
  
  data_split=split(cell_type_split[[cell]], cell_type_split[[cell]]$projid)
  
  #First element
  name=names(data_split)[1]
  data3=data_split[[1]]
  data4=data3[,c(3:2341)]
  data5=colMeans(data4)
  data6=as.data.frame(data5)
  colnames(data6)=name
  main_table=data6
  
  for (i in 2:length(data_split)) {
    name=names(data_split)[i]
    data3=data_split[[i]]
    data4=data3[,c(3:2341)]
    data5=colMeans(data4)
    data6=as.data.frame(data5)
    colnames(data6)=name
    main_table=cbind(main_table,data6)
  }
  
  main_table_new=t(main_table)
  main_table_new=as.data.frame(main_table_new)
  main_table_new2=cbind(projid=rownames(main_table_new),main_table_new)
  
  setwd("F:/PFC_429_final_Spring2021/Vanshika")
  metadata=read.csv("dataset_652_basic_04-23-2020.csv")
  metadata1 <- metadata %>%  
    mutate(projid = as.character(projid))
  
  main_table_anno=left_join(main_table_new2,metadata1,by="projid")
  
  data_table=select(main_table_anno,1:2340,"age_death","pmi","gpath")
  

  
  
  #create data frame to store results
  results <- data.frame()
  
  results_coeff_final <- data.frame()
  results_pval_final <- data.frame()
  
  #loop through module scores and variables
  for(var in names(data_table)[c(2:2340)]){
    results_coeff <- data.frame()
    results_pval <- data.frame()
    for(var2 in names(data_table)[2343]){
      # dynamically generate formula
      fmla <- as.formula(paste0(var,"~", var2,"+","pmi","+","age_death"))
      # fit glm model
      fit <- glm(fmla, data=data_table)
      
      # get coefficents of fit
      cfit <- coef(summary(fit))
      
      # create temporary data frame
      df <- data.frame(module_score = var, variable=var2, coefficient = cfit[var2,'Estimate'],Std.Error=cfit[var2,'Std. Error'],p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
      df_coeff <- data.frame(coefficient = cfit[var2,'Estimate'], stringsAsFactors = F)
      rownames(df_coeff)=var
      colnames(df_coeff)=var2
      df_pval <- data.frame(p.value=cfit[var2,'Pr(>|t|)'], stringsAsFactors = F)
      rownames(df_pval)=var
      colnames(df_pval)=var2
      
      # bind rows of temporary data frame to the results data frame
      results <- rbind(results, df)
      
      #################################################
      if (nrow(results_coeff)>0) {
        results_coeff <- cbind(results_coeff, df_coeff)
      }
      if (nrow(results_coeff)==0) {
        results_coeff <- df_coeff
      }
      
      cell_type_name=names(cell_type_split[cell])
      colnames(results_coeff)=cell_type_name
      #################################################
      
      #################################################
      if (nrow(results_pval)>0) {
        results_pval <- cbind(results_pval, df_pval)
      }
      if (nrow(results_pval)==0) {
        results_pval <- df_pval
      }
      
      cell_type_name=names(cell_type_split[cell])
      colnames(results_pval)=cell_type_name
      #################################################
    }
    #################################################
    if (ncol(results_coeff_final)>0) {
      results_coeff_final <- rbind(results_coeff_final, results_coeff)
    }
    if (ncol(results_coeff_final)==0) {
      results_coeff_final <- results_coeff
    }
    
    #################################################
    
    #################################################
    if (ncol(results_pval_final)>0) {
      results_pval_final <- rbind(results_pval_final, results_pval)
    }
    if (ncol(results_pval_final)==0) {
      results_pval_final <- results_pval
    }
    
    #################################################
  }
  
  # table_coeff_Exc3 (all cell types)
  if(length(table_coeff_Exc3)>0){
    table_coeff_Exc3=cbind(table_coeff_Exc3,results_coeff_final)
  }
  
  if(length(table_coeff_Exc3)==0){
    table_coeff_Exc3=results_coeff_final
  }
  
  
  # table_pval_Exc3 (all cell types)
  if(length(table_pval_Exc3)>0){
    table_pval_Exc3=cbind(table_pval_Exc3,results_pval_final)
  }
  
  if(length(table_pval_Exc3)==0){
    table_pval_Exc3=results_pval_final
  }
}


#################################################################################
#################################################################################

#complex_name=read.csv("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/coreComplexes/coreComplexes.csv")
#complex_name=complex_name[complex_name$Organism %like% "Human",]
#colnames(complex_name)[1]="ComplexID"
#complex_name=complex_name[,1:2]
#complex_name$ComplexID=as.character(complex_name$ComplexID)

# pval table 
table_pval_Exc1=cbind(ComplexID=rownames(table_pval_Exc1),table_pval_Exc1)
table_pval_Exc2=cbind(ComplexID=rownames(table_pval_Exc2),table_pval_Exc2)
table_pval_Exc3=cbind(ComplexID=rownames(table_pval_Exc3),table_pval_Exc3)

table_pval_all=left_join(table_pval_Exc1,table_pval_Exc2, by="ComplexID")
table_pval_all=left_join(table_pval_all,table_pval_Exc3, by="ComplexID")

rownames(table_pval_all)=table_pval_all[,1]
table_pval_all=table_pval_all[,-1]

# coeff table
table_coeff_Exc1=cbind(ComplexID=rownames(table_coeff_Exc1),table_coeff_Exc1)
table_coeff_Exc2=cbind(ComplexID=rownames(table_coeff_Exc2),table_coeff_Exc2)
table_coeff_Exc3=cbind(ComplexID=rownames(table_coeff_Exc3),table_coeff_Exc3)

table_coeff_all=left_join(table_coeff_Exc1,table_coeff_Exc2, by="ComplexID")
table_coeff_all=left_join(table_coeff_all,table_coeff_Exc3, by="ComplexID")

rownames(table_coeff_all)=table_coeff_all[,1]
table_coeff_all=table_coeff_all[,-1]


#BH correction
factor=nrow(table_pval_all)*ncol(table_pval_all)
table_pval2=apply(table_pval_all,2,p.adjust,method="BH", n = factor)

#minus log10(pvalues)
logp_table_pval=-log10(table_pval2)
logp_table_pval[logp_table_pval<0]<-0

#Association score
association_score=logp_table_pval*(table_coeff_all/(abs(table_coeff_all)))


#Output folder
setwd("F:/PFC_429_final_Spring2021/Sudhagar/CORUM_module_score_analysis/coreComplexes/glm/gpath_BH/Exc")

write.csv(table_pval2,file="BH_pval_glm_PFC427_Exc_CORUM_gpath_BH.csv")
write.csv(table_coeff_all,file="Regression_coefficient_glm_PFC427_Exc_CORUM_gpath_BH.csv")
write.csv(logp_table_pval,file="BH_minusLOG10pval_glm_PFC427_Exc_CORUM_gpath_BH.csv")
write.csv(association_score,file="association_score_glm_PFC427_Exc_CORUM_gpath_BH.csv")



#################################################################################
#################################################################################

