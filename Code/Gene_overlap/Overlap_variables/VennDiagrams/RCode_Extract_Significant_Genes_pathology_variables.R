library(stringr)
library(qpcR)

celltype_list=c(
  "Ast CHI3L1",                
  "Ast DPP10",                 
  "Ast GRM3",                  
  "Ast",                       
  "CAMs",                      
  "End",                       
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
  "Fib FLRT2",                 
  "Fib SLC4A4",                
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
  "Mic P2RY12",                
  "Mic TPT1",                  
  "Mic",                       
  "Oli",                       
  "OPC",                       
  "Per",                       
  "SMC",                       
  "T cells")

variable_list=list(
  "amyloid",
  "gpath",
  "nft",
  "tangles",
  "plaq_n",
  "plaq_d")

for (k in 1:length(celltype_list)){
  celltype=celltype_list[k]
  for (j in 1:length(variable_list)){
    variable=variable_list[j]
    working_directory=paste0("D:/PFC_427_WS/Vanshika/muscat/Results/",variable)
    setwd(working_directory)
    filename=paste0(celltype,"_",variable,".csv")
    data=read.csv(filename)
    data_up=subset(data, logFC > 0)
    data_up=subset(data_up, p_adj.loc < 0.05)
    data_down=subset(data, logFC < 0)
    data_down=subset(data_down, p_adj.loc < 0.05)
    
    name=variable
    
    up_genes=data_up[,"gene",drop=FALSE]
    down_genes=data_down[,"gene",drop=FALSE]
    names(up_genes)[names(up_genes) == "gene"] <- name
    names(down_genes)[names(down_genes) == "gene"] <- name
    
    if (j==1){
      table_up=up_genes
      table_down=down_genes
    }
    if (j>1){
      table_up=qpcR:::cbind.na(table_up, up_genes)
      table_down=qpcR:::cbind.na(table_down, down_genes)
    }
  }
  setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap/Overlap_variables/Venn_Diagrams")
  output_file_name_up=paste0(celltype,"_up_genes.csv")
  output_file_name_down=paste0(celltype,"_down_genes.csv")
  write.csv(table_up,file=output_file_name_up)
  write.csv(table_down,file=output_file_name_down)
}

