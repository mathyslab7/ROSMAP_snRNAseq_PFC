library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)

setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores/Intersect_with_consensus_signatures")

data=read.csv("Results_table_scores_complex165_166_5464_Oli_up.csv")
rownames(data)=data$gene
data=data[-c(1:2)]

Cor_matrix=data

#Complex heatmap
#To define clusters
library(circlize)
library(RColorBrewer)
col_fun = colorRamp2(c(-20,-16,-12,-8,-4,0,4,8,12,16,20), rev(brewer.pal(11, "RdBu")))
set.seed(42)
hm=Heatmap(Cor_matrix,col=col_fun,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = FALSE,show_row_names = TRUE,column_names_side = c("top"),show_row_dend=TRUE
        ,row_split=8,
        clustering_distance_rows = "pearson",
        clustering_distance_columns  = "pearson",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#D0E5F0","#85C2EA","#DFF2CF","#B3E9AF","#FDD8D8","#F7B6B6","#FEE4C3","#FFCC99","#E4D7EA","#C1A7DC","#FFFFD6","#E8B69B"),col = 0))),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#D0E5F0","#85C2EA","#DFF2CF","#B3E9AF","#FDD8D8","#F7B6B6","#FEE4C3","#FFCC99","#E4D7EA","#C1A7DC","#FFFFD6","#E8B69B"),col = 0))),
        width = ncol(Cor_matrix)*unit(5, "mm"), 
        height = nrow(Cor_matrix)*unit(0.1, "mm"),
        row_title = NULL,
        column_title = NULL,
        #use_raster=TRUE,
        heatmap_legend_param = list(
          title = "correlation"))

row_order=row_order(hm)
table=data.frame()
for (i in 1:length(row_order)){
  data_1=row_order[[i]]
  data2=as.data.frame(data_1)
  names(data2)[names(data2) == "data_1"] <- "index"
  data2['Cluster'] = i
  if (ncol(table)>0){
    table=rbind(table,data2)
  }
  if (ncol(table)==0){
    table=data2
  } 
}

row_labels=hm@row_names_param[["labels"]]

table2=data.frame()
for (j in 1:length(row_labels)){
  gene=row_labels[[j]]
  data3=as.data.frame(gene)
  data3['index']=j
  if (ncol(table2)>0){
    table2=rbind(table2,data3)
  }
  if (ncol(table2)==0){
    table2=data3
  }
}

table3=left_join(table,table2,by="index")
setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores/Intersect_with_consensus_signatures")
#write.csv(table3,file="Clustering_results_Oli_up_complexs_165_166_5464.csv")

#Select labels and indices (adjust the following line)

#labels=subset(table3, index %in% c(79,84,309,1610))

#labels=read.csv("labels.csv")
#at=labels$index
#genes=labels$gene

#ha = rowAnnotation(foo = anno_mark(at = at, labels = genes))
set.seed(42)

#Original paired color palette
#setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores/Intersect_with_consensus_signatures")
#pdf(file = "CorMatrix_Exc_up.pdf", width = 14, height = 14)
#Heatmap(Cor_matrix,col=col_fun,show_column_names = FALSE,cluster_rows = TRUE,cluster_columns = TRUE,show_row_names = TRUE,column_names_side = c("top"),show_row_dend=TRUE
#        ,row_split=12,column_split=12,
#        clustering_distance_rows = "pearson",
#        clustering_distance_columns  = "pearson",
#        clustering_method_rows = "complete",
#        clustering_method_columns = "complete",
#        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),col = 0))),
#        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),col = 0))),
#        width = ncol(Cor_matrix)*unit(0.1, "mm"), 
#        height = nrow(Cor_matrix)*unit(0.1, "mm"),
#        row_title = NULL,
#        column_title = NULL,
#        use_raster=TRUE,
#        heatmap_legend_param = list(
#          title = "correlation"))+ha

#dev.off()

#Adjusted color palette (paired palette adjusted via https://coolors.co/f9c8c8-efbcd5-d6bddb-8661c1-4b5267) (and rearranged order)
setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores/Intersect_with_consensus_signatures")
pdf(file = "Oli_up_complex165_166_5464_RdBu.pdf", width = 14, height = 14)
Heatmap(Cor_matrix,col=col_fun,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = FALSE,show_row_names = FALSE,column_names_side = c("top"),show_row_dend=TRUE
        ,row_split=8,
        clustering_distance_rows = "pearson",
        clustering_distance_columns  = "pearson",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#D0E5F0","#85C2EA","#DFF2CF","#B3E9AF","#FEE4C3","#FFCC99","#E4D7EA","#C1A7DC","#FFFFD6","#E8B69B"),col = 0))),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#D0E5F0","#85C2EA","#DFF2CF","#B3E9AF","#FEE4C3","#FFCC99","#E4D7EA","#C1A7DC","#FFFFD6","#E8B69B"),col = 0))),
        width = ncol(Cor_matrix)*unit(5, "mm"), 
        height = nrow(Cor_matrix)*unit(0.1, "mm"),
        row_title = NULL,
        column_title = NULL,
        use_raster=FALSE,
        heatmap_legend_param = list(
          title = "signed adj_pval"))

dev.off()



#Code for plotting triangle matrix
#setwd("D:/PFC_427_WS/Vanshika/CORUM_complex/Correlation_results/Scores/Intersect_with_consensus_signatures")
#pdf(file = "GeneOverlap_Oli_up.pdf", width = 14, height = 14)
#Heatmap(Cor_matrix, rect_gp = gpar(type = "none"), column_dend_side = "bottom",
#        col=col_fun,show_column_names = FALSE,cluster_rows = TRUE,cluster_columns = TRUE,show_row_names = TRUE,column_names_side = c("top"),show_row_dend=TRUE
#        ,
#        clustering_distance_rows = "pearson",
#        clustering_distance_columns  = "pearson",
#        clustering_method_rows = "complete",
#        clustering_method_columns = "complete",
#        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),col = 0))),
#        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),col = 0))),
#        width = ncol(Cor_matrix)*unit(0.1, "mm"), 
#        height = nrow(Cor_matrix)*unit(0.1, "mm"),
#        row_title = NULL,
#        column_title = NULL,
#        cell_fun = function(j, i, x, y, w, h, fill) {
#          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#          }
#        })
#dev.off()