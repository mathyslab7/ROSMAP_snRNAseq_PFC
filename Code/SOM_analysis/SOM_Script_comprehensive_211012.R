setwd("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021")
Gene.trait.correlation.matrices.concatenated=read.csv(file="CorMatrix_dplyr_comprehensive.csv",row.names=1)


require(kohonen)
require(RColorBrewer)
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
All_SOM <- som(scale(Gene.trait.correlation.matrices.concatenated), grid = somgrid(26, 26, "hexagonal"), maxNA.fraction=0.5)

save.image("F:/PFC_429_final_Spring2021/Vanshika/SOM_Analysis_Oct2021/SOM_Workspace_211012.RData")