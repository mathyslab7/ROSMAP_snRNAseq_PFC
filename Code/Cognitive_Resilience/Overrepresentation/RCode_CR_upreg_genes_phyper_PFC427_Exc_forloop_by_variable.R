library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(qpcR)

var_list=c(
  "cogn_global_lv_median", "cogng_random_slope_median",
  "gpath_CR_score_median", "plaq_n_CR_score_median", "nft_CR_score_median", "tangles_CR_score_median",
  "gpath_CDR_score_median","plaq_n_CDR_score_median","nft_CDR_score_median","tangles_CDR_score_median"
)


var_table=data.frame()

for (var in 1:length(var_list)) {
  print(var_list[var])
  
  
  set=c("set1","set2","set3")
  
  set_table=data.frame()
  
  for(i in 1:length(set)){
    print(set[i])
    
    setwd("G:/PFC_429_final_Spring2021/Data")
    pbmc = readRDS(paste0("Exc_raw_",set[i],"_final.rds"))
    
    # Add metadata
    data=pbmc@meta.data
    data=data[,"projid",drop=FALSE]
    
    #Add official metadata
    setwd("G:/PFC_429_final_Spring2021/Vanshika")
    metadata=read.csv("dataset_652_basic_04-23-2020.csv")
    metadata = metadata %>% mutate(Apoe_e4 = case_when(
      grepl("22|23|33",apoe_genotype)~"no",
      grepl("24|34|44",apoe_genotype)~"yes"))
    
    # median 
    q50=quantile(na.omit(metadata[,"gpath"]), prob=c(0.5))
    metadata = metadata %>% mutate(gpath_median = 
                                     case_when(gpath <= q50 ~ "low", 
                                               gpath > q50 ~ "high"
                                     ))
    q50=quantile(na.omit(metadata[,"nft"]), prob=c(0.5))
    metadata = metadata %>% mutate(nft_median = 
                                     case_when(nft <= q50 ~ "low", 
                                               nft > q50 ~ "high"
                                     ))
    q50=quantile(na.omit(metadata[,"plaq_n"]), prob=c(0.5))
    metadata = metadata %>% mutate(plaq_n_median = 
                                     case_when(plaq_n <= q50 ~ "low", 
                                               plaq_n > q50 ~ "high"
                                     ))
    q50=quantile(na.omit(metadata[,"tangles"]), prob=c(0.5))
    metadata = metadata %>% mutate(tangles_median = 
                                     case_when(tangles <= q50 ~ "low", 
                                               tangles > q50 ~ "high"
                                     ))
    q50=quantile(na.omit(metadata[,"cogn_global_lv"]), prob=c(0.5))
    metadata = metadata %>% mutate(cogn_global_lv_median = 
                                     case_when(cogn_global_lv <= q50 ~ "low", 
                                               cogn_global_lv > q50 ~ "high"
                                     ))
    q50=quantile(na.omit(metadata[,"cogng_random_slope"]), prob=c(0.5))
    metadata = metadata %>% mutate(cogng_random_slope_median = 
                                     case_when(cogng_random_slope <= q50 ~ "low", 
                                               cogng_random_slope > q50 ~ "high"
                                     ))
    
    # quartiles 
    q25=quantile(na.omit(metadata[,"gpath"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"gpath"]), prob=c(0.75))
    metadata = metadata %>% mutate(gpath_quartile = 
                                     case_when(gpath <= q25 ~ "first_quarter", 
                                               gpath >= q75 ~ "fourth_quarter"
                                     ))
    q25=quantile(na.omit(metadata[,"nft"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"nft"]), prob=c(0.75))
    metadata = metadata %>% mutate(nft_quartile = 
                                     case_when(nft <= q25 ~ "first_quarter", 
                                               nft >= q75 ~ "fourth_quarter"
                                     ))
    q25=quantile(na.omit(metadata[,"plaq_n"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"plaq_n"]), prob=c(0.75))
    metadata = metadata %>% mutate(plaq_n_quartile = 
                                     case_when(plaq_n <= q25 ~ "first_quarter", 
                                               plaq_n >= q75 ~ "fourth_quarter"
                                     ))
    q25=quantile(na.omit(metadata[,"tangles"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"tangles"]), prob=c(0.75))
    metadata = metadata %>% mutate(tangles_quartile = 
                                     case_when(tangles <= q25 ~ "first_quarter", 
                                               tangles >= q75 ~ "fourth_quarter"
                                     ))
    q25=quantile(na.omit(metadata[,"cogn_global_lv"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"cogn_global_lv"]), prob=c(0.75))
    metadata = metadata %>% mutate(cogn_global_lv_quartile = 
                                     case_when(cogn_global_lv <= q25 ~ "first_quarter", 
                                               cogn_global_lv >= q75 ~ "fourth_quarter"
                                     ))
    q25=quantile(na.omit(metadata[,"cogng_random_slope"]), prob=c(0.25))
    q75=quantile(na.omit(metadata[,"cogng_random_slope"]), prob=c(0.75))
    metadata = metadata %>% mutate(cogng_random_slope_quartile = 
                                     case_when(cogng_random_slope <= q25 ~ "first_quarter", 
                                               cogng_random_slope >= q75 ~ "fourth_quarter"
                                     ))
    
    
    
    data2=left_join(data,metadata,by="projid")
    row.names(data2)=row.names(data)
    pbmc=AddMetaData(pbmc,data2)
    
    #Add CDR score variables
    setwd("G:/PFC_429_final_Spring2021/Vanshika")
    CDRscores=read.csv("PFC427_CDR_scores.csv")
    #Add CR score variables
    CRscores=read.csv("CR_scores.csv")
    #Add tangles CR score variable
    tang_CRscores=read.csv("PFC427_tangles_CR_score.csv")
    
    CR_scores=left_join(CDRscores,CRscores, by="projid")
    CR_scores=left_join(CR_scores,tang_CRscores, by="projid")
    CR_scores=na.omit(CR_scores)
    
    
    # median 
    q50=quantile(CR_scores[,"gpath_CR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(gpath_CR_score_median = 
                                       case_when(gpath_CR_score <= q50 ~ "low", 
                                                 gpath_CR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"nft_CR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(nft_CR_score_median = 
                                       case_when(nft_CR_score <= q50 ~ "low", 
                                                 nft_CR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"plaq_n_CR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(plaq_n_CR_score_median = 
                                       case_when(plaq_n_CR_score <= q50 ~ "low", 
                                                 plaq_n_CR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"tangles_CR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(tangles_CR_score_median = 
                                       case_when(tangles_CR_score <= q50 ~ "low", 
                                                 tangles_CR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"gpath_CDR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(gpath_CDR_score_median = 
                                       case_when(gpath_CDR_score <= q50 ~ "low", 
                                                 gpath_CDR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"nft_CDR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(nft_CDR_score_median = 
                                       case_when(nft_CDR_score <= q50 ~ "low", 
                                                 nft_CDR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"plaq_n_CDR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(plaq_n_CDR_score_median = 
                                       case_when(plaq_n_CDR_score <= q50 ~ "low", 
                                                 plaq_n_CDR_score > q50 ~ "high"
                                       ))
    
    q50=quantile(CR_scores[,"tangles_CDR_score"], prob=c(0.5))
    CR_scores = CR_scores %>% mutate(tangles_CDR_score_median = 
                                       case_when(tangles_CDR_score <= q50 ~ "low", 
                                                 tangles_CDR_score > q50 ~ "high"
                                       ))
    
    # quartiles 
    q25=quantile(CR_scores[,"gpath_CR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"gpath_CR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(gpath_CR_score_quartile = 
                                       case_when(gpath_CR_score <= q25 ~ "first_quarter", 
                                                 gpath_CR_score >= q75 ~ "fourth_quarter"
                                       ))
    
    q25=quantile(CR_scores[,"nft_CR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"nft_CR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(nft_CR_score_quartile = 
                                       case_when(nft_CR_score <= q25 ~ "first_quarter", 
                                                 nft_CR_score >= q75 ~ "fourth_quarter"
                                       ))
    q25=quantile(CR_scores[,"plaq_n_CR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"plaq_n_CR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(plaq_n_CR_score_quartile = 
                                       case_when(plaq_n_CR_score <= q25 ~ "first_quarter", 
                                                 plaq_n_CR_score >= q75 ~ "fourth_quarter"
                                       ))
    q25=quantile(CR_scores[,"tangles_CR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"tangles_CR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(tangles_CR_score_quartile = 
                                       case_when(tangles_CR_score <= q25 ~ "first_quarter", 
                                                 tangles_CR_score >= q75 ~ "fourth_quarter"
                                       ))
    
    q25=quantile(CR_scores[,"gpath_CDR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"gpath_CDR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(gpath_CDR_score_quartile = 
                                       case_when(gpath_CDR_score <= q25 ~ "first_quarter", 
                                                 gpath_CDR_score >= q75 ~ "fourth_quarter"
                                       ))
    
    q25=quantile(CR_scores[,"nft_CDR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"nft_CDR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(nft_CDR_score_quartile = 
                                       case_when(nft_CDR_score <= q25 ~ "first_quarter", 
                                                 nft_CDR_score >= q75 ~ "fourth_quarter"
                                       ))
    q25=quantile(CR_scores[,"plaq_n_CDR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"plaq_n_CDR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(plaq_n_CDR_score_quartile = 
                                       case_when(plaq_n_CDR_score <= q25 ~ "first_quarter", 
                                                 plaq_n_CDR_score >= q75 ~ "fourth_quarter"
                                       ))
    q25=quantile(CR_scores[,"tangles_CDR_score"], prob=c(0.25))
    q75=quantile(CR_scores[,"tangles_CDR_score"], prob=c(0.75))
    CR_scores = CR_scores %>% mutate(tangles_CDR_score_quartile = 
                                       case_when(tangles_CDR_score <= q25 ~ "first_quarter", 
                                                 tangles_CDR_score >= q75 ~ "fourth_quarter"
                                       ))
    
    data2=left_join(data,CR_scores,by="projid")
    row.names(data2)=row.names(data)
    pbmc=AddMetaData(pbmc,data2)
    
    
    ################################################################
    
    Idents(pbmc) <- "cell_type_high_resolution"
    Exc_list=sort(levels(pbmc))
    
    Exc_table=data.frame()
    
    for(e in 1:length(Exc_list)){
      print(Exc_list[e])
      
      pbmc_2=subset(pbmc,idents=Exc_list[e])
      DefaultAssay(pbmc_2)<-"RNA"
      
      # Enrichment of CR genes expressing cells - phyper
      
      # phyper - Overrepresenation 
      gene_list=c("HES4",
                  "PDE10A",
                  "RPH3A",
                  "ST6GAL2",
                  "UST"
      )
      
      
      gene_table=data.frame()
      
      for(g in 1:length(gene_list)) {
        print(gene_list[g])
        
        # m: total cells expressing gene in all variable groups (total population)
        counts=pbmc_2@assays[["RNA"]]@counts
        counts=as.data.frame(counts[gene_list[g],drop=FALSE,])
        counts=counts[,colSums(counts)>0]
        counts=as.data.frame(t(counts))
        m=nrow(counts)
        
        
        # m+n: total cells in all variable groups (total population)
        #m_plus_n=nrow(na.omit(pbmc_2@meta.data[,var_list[var],drop=FALSE]))
        m_plus_n=ncol(pbmc_2)
        
        # n: total cells in all variable groups (m+n) - total cells expressing gene in all variable groups (m)
        
        # k: total cells in variable group of interest
        #metadata=pbmc_2@meta.data[pbmc_2@meta.data[,var_list[var]] %like% c("yes|high|fourth_quarter"),]
        #k=nrow(metadata)
        
        metadata=pbmc_2@meta.data[,var_list[var], drop=FALSE]
        k=as.data.frame(table(unlist(metadata)))
        k$Var1=paste0(var_list[var],"_",k$Var1)
        colnames(k)[2] <- "k"
        
        # q: cells expressing gene in variable group 1 (of interest)
        counts_2=cbind(cell_barcodes=rownames(counts),counts)
        metadata=cbind(cell_barcodes=rownames(metadata),metadata)
        metadata=metadata[,c("cell_barcodes",var_list[var])]
        counts_2=left_join(counts_2,metadata,by="cell_barcodes")
        counts_2=counts_2[,-c(1,2)]
        
        q=as.data.frame(table(unlist(counts_2)))
        q$Var1=paste0(var_list[var],"_",q$Var1)
        colnames(q)[2] <- "q"
        if(var_list[var] %like% "median") q=q[order(q[,"Var1"],decreasing = TRUE),]
        
        enrichment=left_join(q,k, by="Var1")
        enrichment=cbind(enrichment,m,m_plus_n)
        
        # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        enrichment$phyper=phyper(enrichment$q, enrichment$m, enrichment$m_plus_n-enrichment$m, enrichment$k, lower.tail = FALSE, log.p = FALSE)
        
        enrichment=enrichment[enrichment$Var1 %like% "median_high",]
        enrichment=enrichment[,"phyper",drop=FALSE]
        row.names(enrichment)=gene_list[g]
        colnames(enrichment)=Exc_list[e]
        
        if(length(gene_table)>0) gene_table=rbind(gene_table,enrichment)
        if(length(gene_table)==0) gene_table=enrichment
        
        
      }
      
      gene_table=cbind(genes=rownames(gene_table),gene_table)
      
      if(length(Exc_table)>0) Exc_table=left_join(Exc_table,gene_table,by="genes")
      if(length(Exc_table)==0) Exc_table=gene_table
      
    }
    
    
    if(length(set_table)>0) set_table=left_join(set_table,Exc_table,by="genes")
    if(length(set_table)==0) set_table=Exc_table
    
  }
  
  
  row.names(set_table)=set_table[,1]
  set_table=set_table[,-1]
  
  
  #BH correction
  factor=nrow(set_table)*ncol(set_table)
  set_table=apply(set_table,2,p.adjust,method="BH", n = factor)
  
  #minus log10(pvalues)
  logp_table_pval=-log10(set_table)
  logp_table_pval[logp_table_pval<0]<-0
  
  assign(var_list[var],logp_table_pval)
  
  setwd("G:/PFC_429_final_Spring2021/Ghada/Cognitive Resilience/phyper/PFC427/Variables")
  write.csv(logp_table_pval,file=paste0("PFC427_Exc_phyper_CR_genes_BH_logp_",var_list[var],"_high.csv"))
  
}



