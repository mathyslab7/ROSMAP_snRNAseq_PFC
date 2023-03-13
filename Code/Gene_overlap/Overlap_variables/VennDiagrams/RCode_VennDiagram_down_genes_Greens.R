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

for (i in 1:length(celltype_list)){
  celltype=celltype_list[i]
  filename=paste0(celltype,"_down_genes.csv")
setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap/Overlap_variables/Venn_Diagrams")
data=read.csv(filename,row.names=1)
data=data[c("gpath","plaq_n","nft","tangles")]

my_list2 <- list()
for(i in 1:ncol(data)) {             # Using for-loop to add columns to list
  my_list2[[i]] <- data[ , i]
  names(my_list2)[i]=colnames(data[i])
}

my_list3=lapply(my_list2, function(x) x[!is.na(x)])

library(ggvenn)
library("ggVennDiagram")
if (length(my_list3[["gpath"]])>0 & length(my_list3[["nft"]])>0 & length(my_list3[["tangles"]])>0 & length(my_list3[["plaq_n"]])>0){
ggVennDiagram(my_list3, label_alpha = 0)

#option1
library(RColorBrewer)
palette_grey <- colorRampPalette(colors = c("grey", "grey"))(6)
scales::show_col(palette_grey)

#option2
library(RColorBrewer)
black_palette = colorRampPalette(c("black", "black"), space = "Lab")
my_black = black_palette(20)

library(ggplot2)
#p <- ggVennDiagram(my_list3, label_alpha = 0,edge_lty = "dashed", edge_size =1)
p <- ggVennDiagram(my_list3, label_alpha = 0,edge_lty = "solid", edge_size =0.5)
# Red Blue
#p + scale_fill_distiller(palette = "RdBu",direction = -1) + scale_color_brewer(palette = "Paired")
#p + scale_fill_distiller(palette = "RdBu",direction = -1) + scale_colour_manual(values=palette_grey)

setwd("D:/PFC_427_WS/Vanshika/muscat/Gene_overlap/Overlap_variables/Venn_Diagrams/Plots")
output_filename=paste0(celltype,"_down_genes.pdf")
pdf(file = output_filename, width = 5, height =5*9/16)
print(p + scale_fill_distiller(palette = "Greens",direction = 1) + scale_colour_manual(values=my_black))
dev.off()
}
}