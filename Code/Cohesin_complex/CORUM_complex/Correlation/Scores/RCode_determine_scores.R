library(dplyr)
setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results")


celltypes=c(
  "Ast",
  "Mic",
  "Oli",
  "Exc L2-3 CBLN2 LINC02306",
  "Exc L3-4 RORB CUX2",
  "Exc L3-5 RORB PLCH1",
  "Exc L4-5 RORB GABRG1",
  "Exc L4-5 RORB IL1RAPL2",
  "Exc L5-6 RORB LINC02196",
  "Exc L5 ET",
  "Exc L5_6 IT Car3",
  "Exc L5_6 NP",
  "Exc L6 CT",
  "Exc L6 THEMIS NFIA",
  "Exc L6b",
  "Exc NRGN",
  "Exc RELN CHD7",
  "Inh ALCAM TRPM3",
  "Inh CUX2 MSR1",
  "Inh ENOX2 SPHKAP",
  "Inh FBN2 EPB41L4A",
  "Inh GPC5 RIT2",
  "Inh L1-2 PAX6 SCGN",
  "Inh L1-6 LAMP5 CA13",
  "Inh L1 PAX6 CA4",
  "Inh L3-5 SST MAFB",
  "Inh L5-6 PVALB STON2",
  "Inh L5-6 SST TH",
  "Inh L6 SST NPY",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh LAMP5 RELN",
  "Inh PTPRK FAM19A1",
  "Inh PVALB CA8 (Chandelier)",
  "Inh PVALB HTR4",
  "Inh PVALB SULF1",
  "Inh RYR3 TSHZ2",
  "Inh SGCD PDE3A",
  "Inh SORCS1 TTN",
  "Inh VIP ABI3BP",
  "Inh VIP CLSTN2",
  "Inh VIP THSD7B",
  "Inh VIP TSHZ2",
  "OPC",
  "End",
  "Per",
  "Fib",
  "SMC")

for (i in 1:length(celltypes)){
  name=celltypes[[i]]
  file1=paste0(name,"_results_table_P.csv")
  file2=paste0(name,"_results_table_r.csv")
  
  setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results")
  Pvals=read.csv(file1,row.names=1)
  rvals=read.csv(file2,row.names=1)
  
  
  X165=Pvals[1,]
  X166=Pvals[2,]
  X5464=Pvals[3,]
  
  X165t=t(X165)
  X166t=t(X166)
  X5464t=t(X5464) 
  
  X165t_adj=p.adjust(X165t[,"X1651"], method = "BH", n = length(X165t[,"X1651"]))
  X165t_adj=as.data.frame(X165t_adj)
  colnames(X165t_adj)[colnames(X165t_adj) == "X165t_adj"] <- "X1651"
  X165_NA=X165t_adj
  X165_NA[X165t_adj==0]<-NA
  min=min(X165_NA,na.rm=TRUE)
  X165t_adj[X165t_adj==0]<-min
  
  
  X166t_adj=p.adjust(X166t[,"X1661"], method = "BH", n = length(X166t[,"X1661"]))
  X166t_adj=as.data.frame(X166t_adj)
  colnames(X166t_adj)[colnames(X166t_adj) == "X166t_adj"] <- "X1661"
  X166_NA=X166t_adj
  X166_NA[X166t_adj==0]<-NA
  min=min(X166_NA,na.rm=TRUE)
  X166t_adj[X166t_adj==0]<-min
  
  X5464t_adj=p.adjust(X5464t[,"X54641"], method = "BH", n = length(X5464t[,"X54641"]))
  X5464t_adj=as.data.frame(X5464t_adj)
  colnames(X5464t_adj)[colnames(X5464t_adj) == "X5464t_adj"] <- "X54641"
  X5464_NA=X5464t_adj
  X5464_NA[X5464t_adj==0]<-NA
  min=min(X5464_NA,na.rm=TRUE)
  X5464t_adj[X5464t_adj==0]<-min
  
  Pvals_adjust=cbind(X165t_adj,X166t_adj,X5464t_adj)
  rvals_adjusted=t(rvals)
  
  Score=(-log10(Pvals_adjust))*(rvals_adjusted/(abs(rvals_adjusted)))
  setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores")
  name=gsub("/", "_", name)
  filename=paste0(name,"_score.csv")
  filename2=paste0(name,"_adj_pval.csv")
  filename3=paste0(name,"_r.csv")
  write.csv(Score,file=filename)
  write.csv(Pvals_adjust,file=filename2)
  write.csv(rvals_adjusted,file=filename3)
}
