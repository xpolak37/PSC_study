a <- metadata_final %>% group_by(Patient) %>%
  group_by(Patient) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  as.data.frame()

nrow(a)

nrow(alpha_ileum_metadata)

a <- spracovanie_tabuliek(filt_colon_asv_tab,filt_colon_taxa_tab)
saveWorkbook(a,paste0("all_levels_colon_filt.xlsx"), overwrite = TRUE)

# supplements ileum
ileum_genus <- aggregate_taxa(filt_ileum_asv_tab,filt_ileum_taxa_tab,taxonomic_level = "Genus")
ileum_genus <- ileum_genus[[1]]

new_df <- data.frame(SeqID=ileum_genus$SeqID)

# norway
ileum_norway_metadata <- filt_ileum_metadata[filt_ileum_metadata$Country !="CZ",]
ileum_genus_norway <- ileum_genus[,c("SeqID",ileum_norway_metadata$SampleID)]

ileum_genus_norway[,-1] <- apply(ileum_genus_norway[,-1],2,function(x) x/sum(x))
all(new_df$SeqID == ileum_genus_norway$SeqID)

## detected
new_df$detected_NOR <- apply(ileum_genus_norway[,-1],1,function(x) any(x>0))

## healthy NOR 
ileum_norway_metadata_healthy <- ileum_norway_metadata[ileum_norway_metadata$Group == "healthy",]
ileum_genus_norway_healthy <- ileum_genus_norway[,c("SeqID",ileum_norway_metadata_healthy$SampleID)]
new_df$detected_NOR_Healthy <- apply(ileum_genus_norway_healthy[,-1],1,function(x) any(x>0))
new_df$NOR_healthy <- t(apply(ileum_genus_norway_healthy[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_healthy_prev <- (rowSums(ileum_genus_norway_healthy[,-1]>0))/(ncol(ileum_genus_norway_healthy[,-1]))

## pre_ltx NOR 
ileum_norway_metadata_preltx <- ileum_norway_metadata[ileum_norway_metadata$Group == "pre_ltx",]
ileum_genus_norway_preltx <- ileum_genus_norway[,c("SeqID",ileum_norway_metadata_preltx$SampleID)]
new_df$detected_NOR_preltx <- apply(ileum_genus_norway_preltx[,-1],1,function(x) any(x>0))
new_df$NOR_preltx <- t(apply(ileum_genus_norway_preltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_preltx_prev <- (rowSums(ileum_genus_norway_preltx[,-1]>0))/(ncol(ileum_genus_norway_preltx[,-1]))

## post_ltx NOR 
ileum_norway_metadata_postltx <- ileum_norway_metadata[ileum_norway_metadata$Group == "post_ltx",]
ileum_genus_norway_postltx <- ileum_genus_norway[,c("SeqID",ileum_norway_metadata_postltx$SampleID)]
new_df$detected_NOR_postltx <- apply(ileum_genus_norway_postltx[,-1],1,function(x) any(x>0))
new_df$NOR_postltx <- t(apply(ileum_genus_norway_postltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_postltx_prev <- (rowSums(ileum_genus_norway_postltx[,-1]>0))/(ncol(ileum_genus_norway_postltx[,-1]))

# CZECH
ileum_czech_metadata <- filt_ileum_metadata[filt_ileum_metadata$Country == "CZ",]
ileum_genus_czech <- ileum_genus[,c("SeqID",ileum_czech_metadata$SampleID)]

ileum_genus_czech[,-1] <- apply(ileum_genus_czech[,-1],2,function(x) x/sum(x))
all(new_df$SeqID == ileum_genus_czech$SeqID)

## detected
new_df$detected_CZ <- apply(ileum_genus_czech[,-1],1,function(x) any(x>0))

## healthy CZ 
ileum_czech_metadata_healthy <- ileum_czech_metadata[ileum_czech_metadata$Group == "healthy",]
ileum_genus_czech_healthy <- ileum_genus_czech[,c("SeqID",ileum_czech_metadata_healthy$SampleID)]
new_df$detected_CZ_Healthy <- apply(ileum_genus_czech_healthy[,-1],1,function(x) any(x>0))
new_df$CZ_healthy <- t(apply(ileum_genus_czech_healthy[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_healthy_prev <- (rowSums(ileum_genus_czech_healthy[,-1]>0))/(ncol(ileum_genus_czech_healthy[,-1]))

## pre_ltx CZ 
ileum_czech_metadata_preltx <- ileum_czech_metadata[ileum_czech_metadata$Group == "pre_ltx",]
ileum_genus_czech_preltx <- ileum_genus_czech[,c("SeqID",ileum_czech_metadata_preltx$SampleID)]
new_df$detected_CZ_preltx <- apply(ileum_genus_czech_preltx[,-1],1,function(x) any(x>0))
new_df$CZ_preltx <- t(apply(ileum_genus_czech_preltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_preltx_prev <- (rowSums(ileum_genus_czech_preltx[,-1]>0))/(ncol(ileum_genus_czech_preltx[,-1]))

## post_ltx CZ 
ileum_czech_metadata_postltx <- ileum_czech_metadata[ileum_czech_metadata$Group == "post_ltx",]
ileum_genus_czech_postltx <- ileum_genus_czech[,c("SeqID",ileum_czech_metadata_postltx$SampleID)]
new_df$detected_CZ_postltx <- apply(ileum_genus_czech_postltx[,-1],1,function(x) any(x>0))
new_df$CZ_postltx <- t(apply(ileum_genus_czech_postltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_postltx_prev <- (rowSums(ileum_genus_czech_postltx[,-1]>0))/(ncol(ileum_genus_czech_postltx[,-1]))

new_df <- new_df[,c("SeqID","detected_NOR","detected_CZ",
                    "detected_NOR_Healthy","NOR_healthy","NOR_healthy_prev",
                    "detected_CZ_Healthy","CZ_healthy","CZ_healthy_prev",
                    "detected_NOR_preltx","NOR_preltx","NOR_preltx_prev",
                    "detected_CZ_preltx","CZ_preltx","CZ_preltx_prev",
                    "detected_NOR_postltx","NOR_postltx","NOR_postltx_prev",
                    "detected_CZ_postltx","CZ_postltx","CZ_postltx_prev")]


write.csv(new_df,"Projects/PSC_study/analysis/results/Q1/supplement3.csv")

# supplements colon
colon_genus <- aggregate_taxa(filt_colon_asv_tab,filt_colon_taxa_tab,taxonomic_level = "Genus")
colon_genus <- colon_genus[[1]]

new_df <- data.frame(SeqID=colon_genus$SeqID)

# norway
colon_norway_metadata <- filt_colon_metadata[filt_colon_metadata$Country !="CZ",]
colon_genus_norway <- colon_genus[,c("SeqID",colon_norway_metadata$SampleID)]

colon_genus_norway[,-1] <- apply(colon_genus_norway[,-1],2,function(x) x/sum(x))
all(new_df$SeqID == colon_genus_norway$SeqID)

## detected
new_df$detected_NOR <- apply(colon_genus_norway[,-1],1,function(x) any(x>0))

## healthy NOR 
colon_norway_metadata_healthy <- colon_norway_metadata[colon_norway_metadata$Group == "healthy",]
colon_genus_norway_healthy <- colon_genus_norway[,c("SeqID",colon_norway_metadata_healthy$SampleID)]
new_df$detected_NOR_Healthy <- apply(colon_genus_norway_healthy[,-1],1,function(x) any(x>0))
new_df$NOR_healthy <- t(apply(colon_genus_norway_healthy[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_healthy_prev <- (rowSums(colon_genus_norway_healthy[,-1]>0))/(ncol(colon_genus_norway_healthy[,-1]))

## pre_ltx NOR 
colon_norway_metadata_preltx <- colon_norway_metadata[colon_norway_metadata$Group == "pre_ltx",]
colon_genus_norway_preltx <- colon_genus_norway[,c("SeqID",colon_norway_metadata_preltx$SampleID)]
new_df$detected_NOR_preltx <- apply(colon_genus_norway_preltx[,-1],1,function(x) any(x>0))
new_df$NOR_preltx <- t(apply(colon_genus_norway_preltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_preltx_prev <- (rowSums(colon_genus_norway_preltx[,-1]>0))/(ncol(colon_genus_norway_preltx[,-1]))

## post_ltx NOR 
colon_norway_metadata_postltx <- colon_norway_metadata[colon_norway_metadata$Group == "post_ltx",]
colon_genus_norway_postltx <- colon_genus_norway[,c("SeqID",colon_norway_metadata_postltx$SampleID)]
new_df$detected_NOR_postltx <- apply(colon_genus_norway_postltx[,-1],1,function(x) any(x>0))
new_df$NOR_postltx <- t(apply(colon_genus_norway_postltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$NOR_postltx_prev <- (rowSums(colon_genus_norway_postltx[,-1]>0))/(ncol(colon_genus_norway_postltx[,-1]))

# CZECH
colon_czech_metadata <- filt_colon_metadata[filt_colon_metadata$Country == "CZ",]
colon_genus_czech <- colon_genus[,c("SeqID",colon_czech_metadata$SampleID)]

colon_genus_czech[,-1] <- apply(colon_genus_czech[,-1],2,function(x) x/sum(x))
all(new_df$SeqID == colon_genus_czech$SeqID)

## detected
new_df$detected_CZ <- apply(colon_genus_czech[,-1],1,function(x) any(x>0))

## healthy CZ 
colon_czech_metadata_healthy <- colon_czech_metadata[colon_czech_metadata$Group == "healthy",]
colon_genus_czech_healthy <- colon_genus_czech[,c("SeqID",colon_czech_metadata_healthy$SampleID)]
new_df$detected_CZ_Healthy <- apply(colon_genus_czech_healthy[,-1],1,function(x) any(x>0))
new_df$CZ_healthy <- t(apply(colon_genus_czech_healthy[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_healthy_prev <- (rowSums(colon_genus_czech_healthy[,-1]>0))/(ncol(colon_genus_czech_healthy[,-1]))

## pre_ltx CZ 
colon_czech_metadata_preltx <- colon_czech_metadata[colon_czech_metadata$Group == "pre_ltx",]
colon_genus_czech_preltx <- colon_genus_czech[,c("SeqID",colon_czech_metadata_preltx$SampleID)]
new_df$detected_CZ_preltx <- apply(colon_genus_czech_preltx[,-1],1,function(x) any(x>0))
new_df$CZ_preltx <- t(apply(colon_genus_czech_preltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_preltx_prev <- (rowSums(colon_genus_czech_preltx[,-1]>0))/(ncol(colon_genus_czech_preltx[,-1]))

## post_ltx CZ 
colon_czech_metadata_postltx <- colon_czech_metadata[colon_czech_metadata$Group == "post_ltx",]
colon_genus_czech_postltx <- colon_genus_czech[,c("SeqID",colon_czech_metadata_postltx$SampleID)]
new_df$detected_CZ_postltx <- apply(colon_genus_czech_postltx[,-1],1,function(x) any(x>0))
new_df$CZ_postltx <- t(apply(colon_genus_czech_postltx[,-1],1,function(x) quantile(x, probs = c(0.25, 0.5, 0.75))))
new_df$CZ_postltx_prev <- (rowSums(colon_genus_czech_postltx[,-1]>0))/(ncol(colon_genus_czech_postltx[,-1]))

new_df <- new_df[,c("SeqID","detected_NOR","detected_CZ",
                    "detected_NOR_Healthy","NOR_healthy","NOR_healthy_prev",
                    "detected_CZ_Healthy","CZ_healthy","CZ_healthy_prev",
                    "detected_NOR_preltx","NOR_preltx","NOR_preltx_prev",
                    "detected_CZ_preltx","CZ_preltx","CZ_preltx_prev",
                    "detected_NOR_postltx","NOR_postltx","NOR_postltx_prev",
                    "detected_CZ_postltx","CZ_postltx","CZ_postltx_prev")]


write.csv(new_df,"Projects/PSC_study/analysis/results/Q1/supplement3.csv")

spracovanie_tabuliek <- function(asv_table, taxa_table){
  res_species <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Species")
  res_genus <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Genus")
  res_family <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Family")
  res_order <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Order")
  res_class <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Class")
  res_phylum <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Phylum")
  res_domain <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Kingdom")
  
  # species
  asv_species <- res_species[[1]]
  taxa_species <- res_species[[2]]
  
  # genus
  asv_genus <- res_genus[[1]]
  taxa_genus <- res_genus[[2]]
  
  # family
  asv_family <- res_family[[1]]
  taxa_family <- res_family[[2]]
  
  # order
  asv_order <- res_order[[1]]
  taxa_order <- res_order[[2]]
  
  # class
  asv_class <- res_class[[1]]
  taxa_class <- res_class[[2]]
  
  # phylum
  asv_phylum <- res_phylum[[1]]
  taxa_phylum <- res_phylum[[2]]
  
  # kingdom
  asv_domain <- res_domain[[1]]
  taxa_domain <- res_domain[[2]]
  
  wb <- createWorkbook(); 
  
  addWorksheet(wb, sheetName ="Kingdom"); 
  writeData(wb, sheet = "Kingdom", asv_domain)
  setColWidths(wb, c("Kingdom"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Phylum"); 
  writeData(wb, sheet = "Phylum", asv_phylum)
  setColWidths(wb, c("Phylum"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Class"); 
  writeData(wb, sheet = "Class", asv_class)
  setColWidths(wb, c("Class"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Order"); 
  writeData(wb, sheet = "Order", asv_order)
  setColWidths(wb, c("Order"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Family"); 
  writeData(wb, sheet = "Family", asv_family)
  setColWidths(wb, c("Family"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Genus"); 
  writeData(wb, sheet = "Genus", asv_genus)
  setColWidths(wb, c("Genus"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Species"); 
  writeData(wb, sheet = "Species", asv_species)
  setColWidths(wb, c("Species"), cols = 1, widths = "auto")
  
  return(wb)
  #saveWorkbook(wb,"IKEM/knihovny/Olivovna/Olivovna.xlsx", overwrite = TRUE)
}

merge_taxa <- function(asv_table,taxa_table,taxonomic_level,names=TRUE){
  taxa_ranks <- colnames(taxa_table)[-which(colnames(taxa_table)=="SeqID")]
  where_level <- which(tolower(taxa_ranks)==tolower(taxonomic_level))
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  taxa_asv_table <- suppressWarnings(taxa_asv_table %>% 
                                       group_by_at(taxa_ranks[1:where_level]) %>%
                                       summarise(across(names(taxa_asv_table)[9:ncol(taxa_asv_table)], sum)))
  
  seq_ids <- apply(taxa_asv_table[,1:where_level],1, function(x){
    a <- paste0(substring(tolower(colnames(taxa_asv_table[,1:where_level])),1,1),"__",x, collapse = ";")
    return(a)
  })
  
  asv_table_sub <- as.data.frame(taxa_asv_table[,(where_level+1):ncol(taxa_asv_table)])
  taxa_table_sub <- as.data.frame(taxa_asv_table[,1:where_level])
  
  asv_table_sub$SeqID <- seq_ids
  taxa_table_sub$SeqID <- seq_ids
  
  asv_table_sub <- asv_table_sub[,c(ncol(asv_table_sub),1:ncol(asv_table_sub)-1)]
  taxa_table_sub <- taxa_table_sub[,c(ncol(taxa_table_sub),1:ncol(taxa_table_sub)-1)]
  return(list(asv_table_sub,taxa_table_sub))
}



############# supplement p-values ------------
######## Q1 -----------

setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q1")
ti <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs pre_ltx")

ti <- ti[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti)[-1] <- paste("healthy vs pre_ltx",colnames(ti)[-1])

ti2 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "pre_ltx vs post_ltx")
ti2 <- ti2[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti2)[-1] <- paste("pre_ltx vs post_ltx",colnames(ti2)[-1])

ti3 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs post_ltx")
ti3 <- ti3[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti3)[-1] <- paste("healthy vs post_ltx",colnames(ti3)[-1])

setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q2")
ti4 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "non-rPSC vs rPSC")
ti4 <- ti4[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti4)[-1] <- paste("non-rPSC vs rPSC",colnames(ti4)[-1])

ti <- merge(ti,ti2,by="SeqID",all = TRUE)
ti <- merge(ti,ti3, by="SeqID", all=TRUE)
ti <- merge(ti,ti4, by="SeqID", all=TRUE)

setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q1")

write.xlsx(ti,"supplement2_pvalues.xlsx")


# COLON
colon <- read.xlsx("univariate_analysis/uni_analysis_wb_colon.xlsx",sheet = "healthy vs pre_ltx")

colon <- colon[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon)[-1] <- paste("healthy vs pre_ltx",colnames(colon)[-1])

colon2 <- read.xlsx("univariate_analysis/uni_analysis_wb_colon.xlsx",sheet = "pre_ltx vs post_ltx")
colon2 <- colon2[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon2)[-1] <- paste("pre_ltx vs post_ltx",colnames(colon2)[-1])

colon3 <- read.xlsx("univariate_analysis/uni_analysis_wb_colon.xlsx",sheet = "healthy vs post_ltx")
colon3 <- colon3[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon3)[-1] <- paste("healthy vs post_ltx",colnames(colon3)[-1])

setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q2")
colon4 <- read.xlsx("univariate_analysis/uni_analysis_wb_colon.xlsx",sheet = "non-rPSC vs rPSC")
colon4 <- colon4[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon4)[-1] <- paste("non-rPSC vs rPSC",colnames(colon4)[-1])


colon <- merge(colon,colon2,by="SeqID",all = TRUE)
colon <- merge(colon,colon3, by="SeqID", all=TRUE)
colon <- merge(colon,colon4, by="SeqID", all=TRUE)
setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q1")

write.xlsx(colon,"supplement3_pvalues.xlsx")

######## Q2 -----------
setwd("C:/Users/povp/Documents/IKEM/Projects/PSC_study/analysis/results/Q2")
ti <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "non-rPSC vs rPSC")

ti <- ti[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti)[-1] <- paste("non-rPSC vs rPSC",colnames(ti)[-1])

ti2 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs non-rPSC")
ti2 <- ti2[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti2)[-1] <- paste("healthy vs non-rPSC",colnames(ti2)[-1])

ti3 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs rPSC")
ti3 <- ti3[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(ti3)[-1] <- paste("healthy vs rPSC",colnames(ti3)[-1])

ti <- merge(ti,ti2,by="SeqID",all = TRUE)
ti <- merge(ti,ti3, by="SeqID", all=TRUE)

write.xlsx(ti,"supplement2_Q2_pvalues.xlsx")


# COLON
colon <- read.xlsx("univariate_analysis/uni_analysis_wb_colon.xlsx",sheet = "non-rPSC vs rPSC")

colon <- colon[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon)[-1] <- paste("non-rPSC vs rPSC",colnames(colon)[-1])

colon2 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs non-rPSC")
colon2 <- colon2[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon2)[-1] <- paste("healthy vs non-rPSC",colnames(colon2)[-1])

colon3 <- read.xlsx("univariate_analysis/uni_analysis_wb_terminal_ileum.xlsx",sheet = "healthy vs rPSC")
colon3 <- colon3[,c("SeqID","log2FoldChange","p_value","padj","final_sig")]
colnames(colon3)[-1] <- paste("healthy vs rPSC",colnames(colon3)[-1])

colon <- merge(colon,colon2,by="SeqID",all = TRUE)
colon <- merge(colon,colon3, by="SeqID", all=TRUE)

write.xlsx(colon,"supplement3_Q2_pvalues.xlsx")
