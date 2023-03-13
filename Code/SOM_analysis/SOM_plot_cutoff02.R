load("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021/SOM_Workspace_211012.RData")
require(kohonen)
require(RColorBrewer)

col=colnames(Gene.trait.correlation.matrices.concatenated)
#col
#setwd("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021")
#write.csv(col,file="column_names_gene_trait_correlation_matrix.csv")

List=as.list(col)

#variable_list

variable_list=list("apoe_genotype",
                   "cogdx",
                   "cogdx_stroke",
                   "dxpark",
                   "cogn_ep_lv",
                   "cogn_po_lv",
                   "cogn_ps_lv",
                   "cogn_se_lv",
                   "cogn_wo_lv",
                   "cognep_random_slope",
                   "cogng_random_slope",
                   "cognpo_random_slope",
                   "cognps_random_slope",
                   "cognse_random_slope",
                   "cognwo_random_slope",
                   "cogn_global_lv",
                   "age_death",
                   "msex",
                   "hypertension_bl",
                   "cancer_bl",
                   "diabetes_sr_rx_bl",
                   "headinjrloc_bl",
                   "heart_bl",
                   "stroke_bl",
                   "bradysc_lv",
                   "gaitsc_lv",
                   "parksc_lv",
                   "braaksc",
                   "ceradsc",
                   "gpath",
                   "gpath_3neocort",
                   "niareagansc",
                   "plaq_d_mf",
                   "pmi",
                   "amyloid",
                   "plaq_d",
                   "plaq_n",
                   "plaq_n_mf",
                   "dlbdx",
                   "nft",
                   "nft_mf",
                   "tangles",
                   "tdp_st4",
                   "arteriol_scler",
                   "caa_4gp",
                   "cvda_4gp2",
                   "ci_num2_gct",
                   "ci_num2_mct",
                   "ci_num2_gtt",
                   "ci_num2_mtt",
                   "gpath_CR_score",
                   "plaq_n_CR_score",
                   "nft_CR_score"
)

#create folders
for (k in seq_along(variable_list)){
  folder<-dir.create(paste0("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021/SOM_plots_cutoff02/",variable_list[k]))
}

for (j in 1:length(variable_list)){
variable=variable_list[[j]]
working_directory=paste0("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021/SOM_plots_cutoff02/",variable)
setwd(working_directory)

library(unikn)
color_Spectral <- c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#FFFFBF","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")
pal_google   <- newpal(color_Spectral)

coolBlueHotRed <- function(n, alpha = 1) {rev(usecol(pal_google, n = n))}

for(i in seq(j, length(List), 53)) {
  name=(List)[[i]]
  pdf(file = paste(name,".pdf",sep = "",collapse = NULL))
  var <- i
  var_unscaled <- aggregate(as.numeric(Gene.trait.correlation.matrices.concatenated[,var]), by=list(All_SOM$unit.classif), FUN=mean, simplify=TRUE,na.rm = TRUE)[,2]
  var_unscaled <- replace(var_unscaled, var_unscaled < -0.2, -0.2)
  var_unscaled <- replace(var_unscaled, var_unscaled > 0.2, 0.2)
  plot(All_SOM, type = "property",ncolors = 40,zlim=c(-0.2,0.2),na.color = "#f4f0ec", property=var_unscaled, main=colnames(getCodes(All_SOM))[var], palette.name=coolBlueHotRed,shape="straight",heatkeywidth = 1,border=FALSE)
  dev.off()
}
}


#Instructions for Combining individual pdf files

#Open Adobe Acrobat DC (icon on Desktop)
#Go to "Tools" and then "Combine files"
#Drag and drop the pdf files you want to combine
#Combine files by clicking on "Combine"
#This will generate a file with 1 pdf per page
#Go to "File" and then "Print"
#Printer: choose "Adobe PDF"
#Under "Page Sizing & Handling" choose "Multiple"
#Under "Pages per sheet" choose "Custom" and the arrangement you want (I think 4 by 5 works pretty well)
#Click "Print", choose the folder and file name, click "Save"

