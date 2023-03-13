library(dplyr)
library(data.table)

setwd("G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/PFC427/muscat/Summary/Gene_Scoring/Upreg")
files=list.files()
files=files[files %like% "Exc"]
files=files[!files %like% "Exc__"]
files=files[!files %like% "Summary"]

all_genes=list()

for(i in 1:length(files)){
  data=read.csv(files[i],row.names = 1)
  data=data[data$count == 16,]
  data=data$name
  name=files[i]
  
  all_genes[[name]]=data
  
}

all_genes=unique(unlist(all_genes))
all_genes=as.data.frame(all_genes)
colnames(all_genes)="name"
all_genes=all_genes[order(all_genes$name),drop=FALSE,]


all_counts=data.frame()
for(i in 1:length(files)){
  data=read.csv(files[i],row.names = 1)
  data=left_join(all_genes,data,by="name")
  data=na.omit(data)
  colnames(data)[2]=gsub("_DEG_upreg_resilience_genes_counts_all_conditions.csv","",files[i])
  
  if(nrow(data)>0){
  if(length(all_counts)>0) all_counts=full_join(all_counts,data,by="name")
  if(length(all_counts)==0) all_counts=data
  }
  
}

rownames(all_counts)=all_counts[,1]
all_counts=all_counts[,-1]
all_counts$count=rowSums(!is.na(all_counts))


write.csv(all_counts,file = "Summary_DEG_upreg_resilience_genes_counts_all_conditions_Exc16.csv")

