suppressMessages(suppressWarnings({
  library(data.table)
  library(cowplot)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(phyloseq)
  library(MicrobiotaProcess)
  library(ggpubr)
  library(ALDEx2)
  library(ggrepel)
  library(ggplotify)
  library(radEmu)
  library(vegan)
  library(reshape2)
  library(pheatmap)
  library(mgcv)
  library(robustlmm)
  library(lmerTest)
  library(emmeans)
  library(magrittr)
  library(openxlsx)
  library(caret)
  library(MicrobiomeStat)
  library(glmnet)
  library(pROC)
  library(purrr)
  library(umap)
  library(Maaslin2)
  library(ggvenn)
  library(ranger)
  library(doParallel)
  library(gbm)
  library(tidyr)
}))
## Basic functions for data processing ----

read_counts <- function(asv_table, text=FALSE, line=5000){
  # This function creates a plot of library size - read count per sample
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # outputs:
  # library size plot
  
  where_seqid <- which(colnames(asv_table)=="SeqID")
  counts <- as.data.frame(colSums(asv_table[,-where_seqid]))
  counts <- counts  %>% setNames(c("Count")) %>% rownames_to_column(var="Sample")
  max_depth <- max(counts$Count)
  
  p_counts <- counts %>%
    ggplot() +
    geom_point( aes(x=reorder(Sample,Count), y=Count), alpha=.8, show.legend = TRUE,colour="#071952") +
    geom_hline(yintercept = line, linetype="dashed", 
               color = "red") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),
      legend.position='none') +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    xlab("Sample") +
    labs(y = "Sequencing depth") + 
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + 
    ylim(0,max_depth)
  #theme_bw() 
  
  if (text) p_counts <- p_counts + 
    geom_text(aes(x=reorder(Sample,Count), y=Count,label = Count), size = 1, hjust = -0.2) 
  return(p_counts)
}

data_check <- function(asv_table,taxa_table){
  # NaN to 0 and unassigned
  asv_table[is.na(asv_table)] <- 0
  taxa_table[is.na(taxa_table)] <- "unassigned"
  taxa_table[taxa_table ==""] <- "unassigned"
  
  rownames(asv_table) <- NULL
  rownames(taxa_table) <- NULL
  
  # check for redundancy
  to_retain <- rowSums(asv_table[,-1]) != 0
  if(FALSE %in% to_retain){
    to_discard <- asv_table$SeqID[!to_retain]
    message((paste("Removing", sum(!to_retain), "ASV(s)")))
    asv_table <- asv_table[-which(asv_table$SeqID %in% to_discard),]
    #taxa_table <- taxa_table[-which(taxa_table$SeqID %in% to_discard),]
  }
  taxa_table %<>% column_to_rownames("SeqID") 
  taxa_table <- taxa_table[asv_table$SeqID,] %>% rownames_to_column("SeqID")
  
  # removing rownames
  row.names(asv_table) <- NULL
  row.names(taxa_table) <- NULL
  
  return(list(asv_table,taxa_table))
}

merging_data <- function(asv_tab_1, asv_tab_2,
                        taxa_tab_1,taxa_tab_2,
                        metadata_1,metadata_2,
                        segment, Q){
  
  if ("segment" %in% colnames(metadata_1)){
    metadata_1 <- metadata_1 %>% dplyr::rename(Matrix=segment)
  }

  # TAB 1
  if (segment=="TI") metadata_1 <- metadata_1[metadata_1$Matrix == segment,]
  else if (segment=="colon") metadata_1 <- metadata_1[metadata_1$Matrix %in% c("Cecum","Rectum","CD","CA","CD","SI"),]
  asv_tab_1 <- asv_tab_1[,c(TRUE,colnames(asv_tab_1)[-1] %in% metadata_1$SampleID)]
  taxa_tab_1 <- taxa_tab_1[taxa_tab_1$SeqID %in% asv_tab_1$SeqID,]
  
  data_checked <- data_check(asv_tab_1,taxa_tab_1)
  asv_tab_1 <- data_checked[[1]]
  taxa_tab_1 <- data_checked[[2]]
  
  if (!is.null(asv_tab_2)){
    # TAB 2
    if (segment=="TI") metadata_2 <- metadata_2[metadata_2$segment == segment,]
    else if (segment=="colon")  metadata_2 <- metadata_2[metadata_2$segment %in% c("CA","CD","SI"),]
    asv_tab_2 <- asv_tab_2[,c(TRUE,colnames(asv_tab_2)[-1] %in% metadata_2$SampleID)]
    taxa_tab_2 <- taxa_tab_2[taxa_tab_2$SeqID %in% taxa_tab_2$SeqID,]
    
    data_checked <- data_check(asv_tab_2,taxa_tab_2)
    asv_tab_2 <- data_checked[[1]]
    taxa_tab_2 <- data_checked[[2]]
    
    # Merging
    merged_asv_tab <- merge(asv_tab_1,asv_tab_2,by="SeqID",all=TRUE)
    merged_taxa_tab <- merging_taxa_tables(taxa_tab_1,taxa_tab_2)
    
    # metadata merge
    merged_metadata <- rbind(metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country),
                             metadata_2 %>% dplyr::select(SampleID,subjectid,Group,segment,Country) %>% dplyr::rename(Patient=subjectid,Matrix=segment) %>%
                               mutate(Patient=paste0("NO_",Patient)))
    row.names(merged_metadata) <- NULL
    
  } else {
    # Merging
    merged_asv_tab <- asv_tab_1
    merged_taxa_tab <- taxa_tab_1
    
    # metadata merge
    if ("Patient" %in% colnames(metadata_1)){
      merged_metadata <- metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country)
    } else {
      merged_metadata <- metadata_1 %>% 
        dplyr::select(SampleID,subjectid,Group,Matrix,Country) %>% 
        dplyr::rename(Patient=subjectid) %>%
        mutate(Patient=paste0("NO_",Patient))
    }
    row.names(merged_metadata) <- NULL
  }
  
  # data check - deleting redundant taxa (all zeros)
  data_checked <- data_check(merged_asv_tab,merged_taxa_tab)
  merged_asv_tab <- data_checked[[1]]
  merged_taxa_tab <- data_checked[[2]]
  
  # keep only Pre_LTx vs Post_LTx vs Healthy
  if (Q == "Q1"){
    merged_metadata <- merged_metadata[merged_metadata$Group %in% c("rPSC","pre_ltx","non-rPSC","healthy"),]
    new_groups <- merged_metadata$Group
    new_groups[merged_metadata$Group %in% c("rPSC","non-rPSC")] <- "post_ltx"
    merged_metadata$Group <- new_groups
  } else if (Q == "Q2"){
    # keep only rPSC vs non-rPSC vs Healthy
    merged_metadata <- merged_metadata[merged_metadata$Group %in% c("rPSC","non-rPSC","healthy"),]
  } else if (Q=="Q3"){
    merged_metadata <- merged_metadata[merged_metadata$Group %in% c("rPSC","pre_ltx"),]
    merged_metadata$Group <- "PSC"
  }
  merged_asv_tab <- merged_asv_tab[,c(TRUE,colnames(merged_asv_tab)[-1] %in% merged_metadata$SampleID)]
  merged_taxa_tab <- merged_taxa_tab[merged_taxa_tab$SeqID %in% merged_asv_tab$SeqID,]
  
  data_checked <- data_check(merged_asv_tab,merged_taxa_tab)
  merged_asv_tab <- data_checked[[1]]
  merged_taxa_tab <- data_checked[[2]]
  row.names(merged_metadata) <- NULL
  
  return(list(merged_asv_tab,merged_taxa_tab, merged_metadata))
}


merging_taxa_tables <- function(taxa_table,new_taxa_table){
  message(("Merging at ASV level"))
  taxa_inconsistent <- merge(taxa_table,new_taxa_table, 
                             by=c("SeqID"), all=TRUE)
  
  # when merging at ASV level, the same sequences can have different taxonomy classification
  # this can happen when they are not classified at once
  # this code finds the inconsistencies in taxonomy and tries to retain the better assignment
  # - the one with less unassigned classification
  
  message(("Finding inconsistencies in taxonomy, trying to keep the ones that have better taxonomy assignment"))
  
  for (i in 1:nrow(taxa_inconsistent)){
    taxonomy_old <- taxa_inconsistent[i,2:8] 
    taxonomy_new <- taxa_inconsistent[i,9:15] 
    #taxonomy_old <- taxa_inconsistent[i,2:7] 
    #taxonomy_new <- taxa_inconsistent[i,8:13] 
    if ((sum(is.na(taxonomy_old))!=7) &(sum(is.na(taxonomy_new))!=7) ){
      if (FALSE %in% (taxonomy_old == taxonomy_new)){
        old_sum <- sum(taxonomy_old == "unassigned")
        new_sum <- sum(taxonomy_new == "unassigned")
        if (old_sum > new_sum){ 
          my_taxonomy <- taxonomy_new
        } else my_taxonomy <- taxonomy_old
        
        # assigning same taxonomy
        new_taxa_table[new_taxa_table$SeqID==(taxa_inconsistent$SeqID[i]),-1] <- my_taxonomy
        taxa_table[taxa_table$SeqID==(taxa_inconsistent$SeqID[i]),-1] <- my_taxonomy
      }
    }
  }
  
  merged_df <- rbind(taxa_table,new_taxa_table)
  
  # keep unique ASVs
  taxa_table <- merged_df[!duplicated(merged_df$SeqID),]
  rownames(taxa_table) <- NULL
  return(taxa_table)
}

aggregate_taxa <- function(asv_table,taxa_table,taxonomic_level,names=TRUE){
  # This function aggregates abundances at different taxonomic levels.
  # It is used in every function with parameter 'taxonomic_level', where
  # various taxonomic levels can be considered.
  # inputs:
  # asv_table - ASV table with 'SeqID'
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # taxonomic_level - aggregate to which level
  # names=TRUE - boolean indicating if 'SeqID' should be returned, several functions
  # do not require the IDs, so this column can be discarded
  # outputs:
  # list with new asv table and taxa table
  # available taxa ranks
  capitalize_first <- function(string) {
    paste0(toupper(substr(string, 1, 1)), tolower(substr(string, 2, nchar(string))))
  }
  taxonomic_level <- capitalize_first(taxonomic_level)
  taxa_ranks <- colnames(taxa_table)[-which(colnames(taxa_table)=="SeqID")]
  where_level <- which(taxa_ranks==taxonomic_level)
  
  # merge asv and taxa table, this table is than used for aggregating
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  taxa_asv_table <- suppressMessages(suppressWarnings(taxa_asv_table %>% 
                                                        group_by_at(taxa_ranks[1:where_level]) %>%
                                                        summarise(across(names(taxa_asv_table)[9:ncol(taxa_asv_table)], sum))))
  
  # splitting the taxa_asv into two dataframes  
  asv_table_sub <- as.data.frame(taxa_asv_table[,(where_level+1):ncol(taxa_asv_table)])
  taxa_table_sub <- as.data.frame(taxa_asv_table[,1:where_level])
  
  # for some purposes, the names can be added (e.g. for visualization), 
  # otherwise names are not provided 
  if (names==TRUE){
    if (taxonomic_level == "Species") asv_table_sub$SeqID <- paste(taxa_table_sub[,where_level-1], taxa_table_sub[,where_level])
    else if (taxonomic_level == "Genus") {
      # genus assigning
      genus <- taxa_table_sub[,where_level]
      
      # family assigning
      where_unassigned <- grep("unassigned",genus)
      where_uncultured <- grep("uncultured",genus)
      where_unassigned <- c(where_unassigned,where_uncultured)
      genus[where_unassigned] <- paste0("f__",taxa_table_sub[where_unassigned,where_level-1],";g__",genus[where_unassigned])
      
      # order assigning
      where_f_unassigned <-  grep("f__unassigned",genus)
      where_f_uncultured <- grep("f__uncultured",genus)
      where_f_unassigned <- c(where_f_unassigned,where_f_uncultured)
      genus[where_f_unassigned] <- paste0("o__",taxa_table_sub[where_f_unassigned,where_level-2],";",genus[where_f_unassigned])
      
      duplicated <- genus[duplicated(genus)]
      
      for (duplic in unique(duplicated)){
        where_duplic <- which(genus==duplic)
        genus[where_duplic] <- paste(duplic,1:length(where_duplic))
      }
      asv_table_sub$SeqID <- genus
    } else asv_table_sub$SeqID <- taxa_table_sub[,where_level]
    
    #duplicated <- asv_table_sub$SeqID[duplicated(asv_table_sub$SeqID)]
    #asv_table_sub$SeqID[duplicated(asv_table_sub$SeqID)] <- paste(duplicated,1:length(duplicated))
    taxa_table_sub$SeqID <- asv_table_sub$SeqID
    
    asv_table_sub <- asv_table_sub[,c("SeqID",colnames(asv_table_sub)[-dim(asv_table_sub)[2]])]
    taxa_table_sub <- taxa_table_sub[,c("SeqID",colnames(taxa_table_sub)[-dim(taxa_table_sub)[2]])]
  }
  return(list(asv_table_sub,taxa_table_sub))
}

aggregate_samples <- function(asv_table,taxa_table,metadata,variable){
  # This function aggregates abundances at different taxonomic levels.
  # It is used in every function with parameter 'taxonomic_level', where
  # various taxonomic levels can be considered.
  # inputs:
  # asv_table - ASV table with 'SeqID'
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # taxonomic_level - aggregate to which level
  # names=TRUE - boolean indicating if 'SeqID' should be returned, several functions
  # do not require the IDs, so this column can be discarded
  # outputs:
  # list with new asv table and taxa table
  # available taxa ranks

  # merge asv and taxa table, this table is than used for aggregating
  group_asv_table <- data.frame(SeqID=taxa_table$SeqID)
  unique_values_variable <- unique(metadata[,variable])
  for (value in unique_values_variable){
      samples <- metadata$SampleID[metadata[,variable]==value]
      group_sums <- rowSums(asv_table[,samples])
      group_asv_table[,value] <- group_sums
  }
  
  return(group_asv_table)
}


construct_phyloseq <- function(asv_table, taxa_table, metadata){
  # This function constructs phyloseq object and is mainly used in functions
  # with methods that require phyloseq object 
  
  otu_mat <- asv_table
  tax_mat <- taxa_table
  samples_df <- metadata
  
  # removing rownames
  rownames(otu_mat) <- NULL
  rownames(tax_mat) <- NULL
  rownames(samples_df) <- NULL
  
  # set SeqID as rownames
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("SeqID") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("SeqID")
  
  # if more than one samples is provided
  if (nrow(samples_df)>1){
    samples_df <- as.data.frame(samples_df) %>% 
      tibble::column_to_rownames("SampleID") 
  } else{
    # if only one sample - rownames have to be set manually
    rownames(samples_df) <- samples_df$SampleID
    samples_df <- samples_df[, -which(names(samples_df) == "SampleID")]
  }
  
  # construct matrix
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
  # preprocessing for phyloseq object
  OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(tax_mat) 
  samples = sample_data(samples_df)
  
  # resulting object
  object <- phyloseq(OTU, TAX, samples)
  return(object)
}

create_asv_taxa_table <- function(asv_table, taxa_table){
  # creating asv+taxa
  asvs <- asv_table$SeqID
  taxa_ranks <- colnames(taxa_table)
  where_level <- which(tolower(taxa_ranks)=="species")
  if (length(where_level)==0) where_level <- which(tolower(taxa_ranks)=="genus")
  if (length(where_level)==0) where_level <- which(tolower(taxa_ranks)=="phylum")
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  
  seq_ids <- apply(taxa_asv_table[,2:where_level],1, function(x){
    a <- paste0(substring(tolower(colnames(taxa_asv_table[,2:where_level])),1,1),"__",x, collapse = ";")
    return(a)
  })
  
  taxa_asv_table$Taxonomy <- seq_ids
  taxa_asv_table %<>% column_to_rownames("SeqID")
  taxa_asv_table <- taxa_asv_table[asvs,]
  taxa_asv_table <- taxa_asv_table[,-which(colnames(taxa_asv_table) %in% taxa_ranks[2:8])]
  taxa_asv_table <- taxa_asv_table[,c(ncol(taxa_asv_table),1:(ncol(taxa_asv_table)-1))]
  colnames(taxa_asv_table) <- c("SeqID", colnames(taxa_asv_table)[-1])
  return(taxa_asv_table)
}

binomial_prep <- function(asv_table,taxa_table,metadata,group, patient=FALSE,
                          usage="linDA"){
  
  # prepare the data
  uni_data <- asv_table %>% column_to_rownames("SeqID")
  metadata %<>% column_to_rownames("SampleID")
  uni_data <- uni_data[,rownames(metadata)[metadata$Group %in% group]]
  
  data_checked <- data_check(
    uni_data %>% rownames_to_column("SeqID"),
    taxa_table)
  
  uni_data <- data_checked[[1]] %>% column_to_rownames("SeqID")
  uni_tax <- data_checked[[2]]
  uni_metadata <- metadata[colnames(uni_data),]
  
  filt_data <- filtering_steps(uni_data %>% rownames_to_column("SeqID"),uni_tax,
                               uni_metadata %>% rownames_to_column("SampleID"),
                               seq_depth_threshold=10000)
  
  uni_data <- filt_data[[1]]  %>% column_to_rownames("SeqID")
  uni_tax <- filt_data[[2]]
  uni_metadata <- filt_data[[3]] %>% column_to_rownames("SampleID")
  
  if (grepl("ml_",usage)){
    if (usage=="ml_clr"){
      uni_data <- vegan::decostand(uni_data,method = "clr", MARGIN = 2,pseudocount=0.5) %>% 
        as.matrix()
    } else if(usage=="ml_ra"){
        uni_data <- as.data.frame(apply(uni_data, 2, function(x) x / sum(x)))
    }
      uni_data <- uni_data %>% t() %>% as.data.frame()
      uni_data$Group <- uni_metadata[rownames(uni_data),"Group"]
      uni_data$Country <- uni_metadata[rownames(uni_data),"Country"]
      if (patient) uni_data$Patient <- uni_metadata[rownames(uni_data),"Patient"]
      
      uni_data %<>%
        mutate(Group = if_else(Group == group[1], 1, 0))
      return(uni_data)
  } else if (usage=="linDA"){
    uni_metadata$Group <- factor(uni_metadata$Group)
    uni_metadata$Group <- relevel(uni_metadata$Group,group[1])
    return(list(uni_data,uni_tax,uni_metadata))
  } else message("Invalid usage - use one of ml_* or linDA")
}


binomial_prep_psc_effect <- function(asv_table,taxa_table,metadata,df_effect, patient=FALSE,
                          usage="linDA"){
  
  if (TRUE %in% grepl("rPSC",unique(metadata$Group))){
    group <- c("rPSC","non-rPSC")
  } else group <- c("pre_ltx","post_ltx")
  
  # prepare the data
  uni_data <- asv_table %>% column_to_rownames("SeqID")
  metadata %<>% column_to_rownames("SampleID")
  uni_data <- uni_data[,rownames(metadata)[metadata$Group %in% group]]
  
  data_checked <- data_check(
    uni_data %>% rownames_to_column("SeqID"),
    taxa_table)
  
  uni_data <- data_checked[[1]] %>% column_to_rownames("SeqID")
  uni_tax <- data_checked[[2]]
  uni_metadata <- metadata[colnames(uni_data),]
  
  filt_data <- seq_depth_filtering(uni_data %>% rownames_to_column("SeqID"),uni_tax,
                               uni_metadata %>% rownames_to_column("SampleID"),
                               seq_depth_threshold=10000)
  
  uni_data <- filt_data[[1]]  %>% column_to_rownames("SeqID")
  uni_tax <- filt_data[[2]]
  uni_metadata <- filt_data[[3]] %>% column_to_rownames("SampleID")
  
  # only selected TAXA
  wanted_asvs <- df_effect$SeqID
  uni_data <- uni_data[wanted_asvs,]
  
  data_checked <- data_check(
    uni_data %>% rownames_to_column("SeqID"),
    uni_tax)
  
  uni_data <- data_checked[[1]] %>% column_to_rownames("SeqID")
  uni_tax <- data_checked[[2]]
  uni_metadata <- metadata[colnames(uni_data),]
  
  if (grepl("ml_",usage)){
    if (usage=="ml_clr"){
      uni_data <- vegan::decostand(uni_data,method = "clr", MARGIN = 2,pseudocount=0.5) %>% 
        as.matrix()
    } else if(usage=="ml_ra"){
      uni_data <- as.data.frame(apply(uni_data, 2, function(x) x / sum(x)))
      uni_data[is.na(uni_data)] <- 0
    }
    
    uni_data <- uni_data %>% t() %>% as.data.frame()
    uni_data$Group <- uni_metadata[rownames(uni_data),"Group"]
    uni_data$Country <- uni_metadata[rownames(uni_data),"Country"]
    if (patient) uni_data$Patient <- uni_metadata[rownames(uni_data),"Patient"]
    
    uni_data %<>%
      mutate(Group = if_else(Group == group[1], 1, 0))
    return(uni_data)
  } else if (usage=="linDA"){
    uni_metadata$Group <- factor(uni_metadata$Group)
    uni_metadata$Group <- relevel(uni_metadata$Group,group[1])
    return(list(uni_data,uni_tax,uni_metadata))
  } else message("Invalid usage - use one of glmnet or linDA")
  
}

univariate_statistics <- function(list_intersections,psc_effect,
                                  genus_asv_taxa_tab,segment="terminal_ileum"){
  univar_df <- data.frame()
  wb = createWorkbook()
  
  # group1 vs group2
  groups_names <- c()
  for (name in (grep(paste(segment,"ASV"),names(list_intersections),value=TRUE))){
    groups_names <- c(groups_names,gsub("_ltx","",gsub("healthy","H",gsub(segment,"",gsub(" ASV","", name)))))
    ASV_df <- list_intersections[[name]]
    genus_df <- list_intersections[[gsub("ASV","Genus",name)]]
    
    ASV <- ASV_df$Taxonomy
    ASV <- substring(ASV, 1, nchar(ASV)-14)
    
    genus <- genus_df$ASV
    genus <- genus_asv_taxa_tab[genus,"SeqID"]
    
    ASVs_in_genera <- ASV %in% genus
    ASVs_in_genera <- ASV_df[ASVs_in_genera,]
    
    unique_ASVs_in_genera <- unique(ASVs_in_genera$Taxonomy)
    unique_ASVs_in_genera <- genus_asv_taxa_tab[unique_ASVs_in_genera,"SeqID"]
    
    ASVs_not_in_genera <- !(ASV %in% genus)
    ASVs_not_in_genera <- ASV_df[ASVs_not_in_genera,]
    
    genera_in_ASVs <- genus %in% unique(ASV)
    genera_in_ASVs <- genus_df[genera_in_ASVs,]
    
    genera_not_in_ASVs <- !(genus %in% unique(ASV))
    genera_not_in_ASVs <- genus_df[genera_not_in_ASVs,]
    
    univar_df <- rbind(univar_df,t(data.frame(
      c(length(ASV),
        length(unique(ASV)), 
        length(genus)
        #nrow(ASVs_in_genera),
        #nrow(ASVs_not_in_genera),
        #length(unique_ASVs_in_genera),
        #nrow(genera_in_ASVs),
        #nrow(genera_not_in_ASVs)
        ))))
    
    new_name <- groups_names[length(groups_names)]
    addWorksheet(wb, sheetName = paste(new_name,"ASVs"))
    writeData(wb, sheet = paste(new_name,"ASVs"), ASV_df, rowNames=FALSE)
    
    addWorksheet(wb, sheetName = paste(new_name,"Genera"))
    writeData(wb, sheet = paste(new_name,"Genera"), genus_df, rowNames=FALSE)
    
    #addWorksheet(wb, sheetName = paste(new_name,"ASVs in G"))
    #writeData(wb, sheet = paste(new_name,"ASVs in G"), ASVs_in_genera, rowNames=FALSE)
    
    #addWorksheet(wb, sheetName = paste(new_name,"ASVs !in G"))
    #writeData(wb, sheet = paste(new_name,"ASVs !in G"), ASVs_not_in_genera, rowNames=FALSE)
    
    #addWorksheet(wb, sheetName = paste(new_name,"AUG in ASVs"))
    #writeData(wb, sheet = paste(new_name,"AUG in ASVs"), unique_ASVs_in_genera, rowNames=FALSE)
    
    #addWorksheet(wb, sheetName = paste(new_name,"G in ASVs"))
    #writeData(wb, sheet = paste(new_name,"G in ASVs"), genera_in_ASVs, rowNames=FALSE)
    
    #addWorksheet(wb, sheetName = paste(new_name,"G !in ASVs"))
    #writeData(wb, sheet = paste(new_name,"G !in ASVs"), genera_in_ASVs, rowNames=FALSE)
    }
  
  # PSC effect/rPSC effect
  groups_names <- c(groups_names,"PSC effect")
  ASV_df <- psc_effect[[paste(segment,"ASV")]]
  genus_df <- psc_effect[[paste(segment,"Genus")]]
  
  ASV <- ASV_df$Taxonomy
  ASV <- substring(ASV, 1, nchar(ASV)-14)
  
  genus <- genus_df$ASV
  genus <- ileum_genus_asv_taxa_tab[genus,"SeqID"]
  
  ASVs_in_genera <- ASV %in% genus
  ASVs_in_genera <- ASV_df[ASVs_in_genera,]
  
  unique_ASVs_in_genera <- unique(ASVs_in_genera$Taxonomy)
    
  ASVs_not_in_genera <- !(ASV %in% genus)
  ASVs_not_in_genera <- ASV_df[ASVs_not_in_genera,]
  
  genera_in_ASVs <- genus %in% unique(ASV)
  genera_in_ASVs <- genus_df[genera_in_ASVs,]
  
  genera_not_in_ASVs <- !(genus %in% unique(ASV))
  genera_not_in_ASVs <- genus_df[genera_not_in_ASVs,]
  
  univar_df <- rbind(univar_df,t(data.frame(
    c(length(ASV),
      length(unique(ASV)), 
      length(genus)
      #nrow(ASVs_in_genera),
      #nrow(ASVs_not_in_genera),
      #length(unique_ASVs_in_genera),
      #nrow(genera_in_ASVs),
      #nrow(genera_not_in_ASVs)
    ))))
  
  rownames(univar_df) <-groups_names
  colnames(univar_df) <- c("ASVs","Unique genera", "Genera")#, 
                           #"ASVs in genera", "ASVs NOT in genera", "Unique ASVs in genera",
                           #"Genera in ASVs",  "Genera NOT in ASVs")
  
  new_name <- groups_names[length(groups_names)]
  addWorksheet(wb, sheetName = paste(new_name,"ASVs"))
  writeData(wb, sheet = paste(new_name,"ASVs"), ASV_df, rowNames=FALSE)
  
  addWorksheet(wb, sheetName = paste(new_name,"Genera"))
  writeData(wb, sheet = paste(new_name,"Genera"), genus_df, rowNames=FALSE)
  
  #addWorksheet(wb, sheetName = paste(new_name,"ASVs in G"))
  #writeData(wb, sheet = paste(new_name,"ASVs in G"), ASVs_in_genera, rowNames=FALSE)
  
  #addWorksheet(wb, sheetName = paste(new_name,"ASVs !in G"))
  #writeData(wb, sheet = paste(new_name,"ASVs !in G"), ASVs_not_in_genera, rowNames=FALSE)
  
  #addWorksheet(wb, sheetName = paste(new_name,"AUG in ASVs"))
  #writeData(wb, sheet = paste(new_name,"AUG in ASVs"), unique_ASVs_in_genera, rowNames=FALSE)
  
  #addWorksheet(wb, sheetName = paste(new_name,"G in ASVs"))
  #writeData(wb, sheet = paste(new_name,"G in ASVs"), genera_in_ASVs, rowNames=FALSE)
  
  #addWorksheet(wb, sheetName = paste(new_name,"G !in ASVs"))
  #writeData(wb, sheet = paste(new_name,"G !in ASVs"), genera_in_ASVs, rowNames=FALSE)
  
  return(list(univar_df, wb))
}

mock_zymo_genus_merging <- function(asv_table,taxa_table,zymo_asv_table, zymo_taxa){
  # merging mock community samples with reference abundance
  taxa_table[is.na(taxa_table)] <- "unassigned"
  asv_table[is.na(asv_table)] <- 0
  taxa_reads_table <- create_asv_taxa_table(asv_table,taxa_table)
  genus_data <- aggregate_taxa(asv_table,taxa_table, taxonomic_level = "Genus")
  genus_asv_table <- genus_data[[1]]
  genus_taxa_table <- genus_data[[2]]
  
  genus_asv_table_norm <- as.data.frame(apply(genus_asv_table[,-1],2,function(x) x/sum(x)))
  genus_asv_table_norm$SeqID <- genus_asv_table$SeqID
  
  colnames(zymo_asv_table)[1] <- "SeqID"
  colnames(zymo_taxa)[1] <- "SeqID"
  
  genus_data_zymo <- aggregate_taxa(zymo_asv_table,zymo_taxa, taxonomic_level = "Genus")
  genus_zymo_asv_table <- genus_data_zymo[[1]][-9,]
  genus_zymo_taxa_table <- genus_data_zymo[[2]][-9,]
  
  merged_data <- merge(genus_zymo_asv_table,genus_asv_table_norm, by="SeqID", all=TRUE)
  merged_data[is.na(merged_data)] <- 0
  colnames(merged_data)[2] <- "ZYMO REFERENCE"
  return(merged_data)
}

group_intersection <- function(group, list_intersections, list_venns,
                                 linda.output1, fit_data,
                                 raw_linda_results, segment,level){
  
  list_core <- list()
  linda_df <- linda.output1[[paste0(group[1]," vs Group",group[2])]]
  list_core[["linDA"]] <- rownames(linda_df[linda_df$padj <= 0.05,])
  
  maaslin_df <- fit_data$results[fit_data$results$metadata=="Group",]
  list_core[["MaAsLin2"]] <- maaslin_df[maaslin_df$qval <= 0.05,"feature"]
  
  diff <- intersect(list_core[[1]],list_core[[2]])
  orig <- raw_linda_results[[segment]][[paste0(group[1]," vs Group",group[2])]]
  diff <- orig[orig$SeqID %in% diff,]
  
  #if (segment == "terminal_ileum") segment <- "Ileum"
  list_intersections[[paste(segment,level,paste(group, collapse = " vs "))]] <- diff
  
  # venn diagram
  venn <- ggvenn(list_core, fill_color = c("blue", "red")) + 
    ggtitle(paste(segment,level,paste(group, collapse = " vs "))) + 
    guides(fill = "none")
  
  # save the results
  list_venns[[paste(segment,level,paste(group, collapse = " vs "))]] <- venn
  
  # show the results
  return(list(list_intersections,list_venns,venn))
}

country_union <- function(group,linda.output1, fit_data, segment, level){
  list_core <- list()
  linda_df <- linda.output1[[paste(group[1],",",group[2],"- CZ vs NO")]]
  list_core[["linDA"]] <- rownames(linda_df[linda_df$padj <= 0.05,])
  
  maaslin_df <- fit_data$results[fit_data$results$metadata=="Group",]
  list_core[["MaAsLin2"]] <- maaslin_df[maaslin_df$qval <= 0.05,"feature"]
  
  union <- c(list_core[[1]],list_core[[2]])
  union <- union[!duplicated(union)]
  
  #if (segment == "terminal_ileum") segment <- "Ileum"
  # save the results
  list_country_union[[paste(segment,level,paste(group, collapse = " vs "))]] <- union
  return(list_country_union)
}

country_interaction <- function(group,linda.output1,list_intersections, 
                                uni_data, uni_metadata,
                                segment, level){
  interaction_significant_df <- linda.output1[[paste0(group[1]," vs Group",group[2],":CountryNO")]]
  interaction_significant_df <- interaction_significant_df[interaction_significant_df$padj<0.05,]
  
  group_significant_df <- list_intersections[[paste(segment,level, group[1],"vs",group[2])]]
  
  intersection_significant <- intersect(rownames(interaction_significant_df),group_significant_df$ASV)
  
  linda_czech_interaction_significant <- NA
  linda_no_interaction_significant <- NA
  
  if (length(intersection_significant)>0){
    # Subset the data by country
    uni_metadata_czech <- subset(uni_metadata, Country == "CZ")
    uni_metadata_no <- subset(uni_metadata, Country == "NO")
    
    uni_data_czech <- uni_data[,rownames(uni_metadata_czech)]
    uni_data_no <- uni_data[,rownames(uni_metadata_no)]
    
    # Run linDA for each subset
    if (segment=="terminal_ileum"){
      linda_czech <- linda(uni_data_czech , uni_metadata_czech, formula = '~ Group')
      linda_no <- linda(uni_data_no, uni_metadata_no, formula = '~ Group')
    } else if (segment=="colon"){
      linda_czech <- linda(uni_data_czech , uni_metadata_czech, formula = '~ Group + (1|Patient)')
      linda_no <- linda(uni_data_no, uni_metadata_no, formula = '~ Group + (1|Patient)')
      
    } else cat("Problem with segment!\n")
    
    linda_czech_interaction_significant <- linda_czech$output[[1]][intersection_significant,]
    linda_no_interaction_significant <- linda_no$output[[1]][intersection_significant,]
  }
  
  return(list(intersection_significant,linda_czech_interaction_significant,linda_no_interaction_significant))
}

removing_interaction_problems <- function(group,
                                          list_interaction_significant,
                                          list_intersections,
                                          segment, level){
  
  to_investigate <- list_interaction_significant[[1]]
  to_remove <- c()
  
  if (length(to_investigate)>0){
    for (taxon in to_investigate){
      cz <- list_interaction_significant[[2]][taxon,c("log2FoldChange","reject")]
      no <- list_interaction_significant[[3]][taxon,c("log2FoldChange","reject")]
      if (! (((cz[,1]>0)==(no[,1]>0)) & cz[,2] & no[,2]) ) to_remove <- c(to_remove,taxon)
    }
    
    if (length(to_remove)>0){
      df_to_change <- list_intersections[[paste(segment,level, group[1],"vs",group[2])]]
      df_to_change <- df_to_change[!(rownames(df_to_change)%in%to_remove),]
      list_intersections[[paste(segment,level, group[1],"vs",group[2])]] <- df_to_change
    }
    
    list_intersections <- lapply(list_intersections, function (x) remove_rownames(x))
  }
  rownames(list_intersections) <- NULL
  return(list_intersections)
}

dysbiosis_index_calculation <- function(my_table, metadata_table, 
                                        psc_increased,
                                        psc_decreased,
                                        name){
  my_table <- my_table %>% column_to_rownames("SeqID")
  
  dysbiosis_data <- data.frame()
  for (i in 1:ncol(my_table)){
    SampleID <- colnames(my_table)[i]
    PatientID <- metadata_table[metadata_table$SampleID==SampleID,"Patient"]
    abundances <- my_table[,i]/sum(my_table[,i])
    names(abundances) <- rownames(my_table)
    abundances_psc_increased <- sum(abundances[psc_increased])
    abundances_psc_decreased <- sum(abundances[psc_decreased])
    
    # !!!!!!!!!!! PSEUDOCOUNT - prediskutovat
    if (abundances_psc_increased==0) abundances_psc_increased <- 1e-20
    if (abundances_psc_decreased==0) abundances_psc_decreased <- 1e-20
    dys_index <- log(abundances_psc_increased/abundances_psc_decreased)
    dysbiosis_data <- rbind(dysbiosis_data,data.frame(SampleID,PatientID,dys_index))
  }
  colnames(dysbiosis_data) <- c("SampleID","PatientID",name)
  return(dysbiosis_data)
}

dysbiosis_index_calculation_clr <- function(my_table, metadata_table, 
                                        psc_increased,
                                        psc_decreased,
                                        name){
  my_table <- my_table %>% column_to_rownames("SeqID")
  data_clr <- vegan::decostand(my_table,method = "clr", MARGIN = 2,pseudocount=0.5) %>% as.matrix()
  
  dysbiosis_data <- data.frame()
  for (i in 1:ncol(my_table)){
    SampleID <- colnames(my_table)[i]
    PatientID <- metadata_table[metadata_table$SampleID==SampleID,"Patient"]
    abundances <- data_clr[,i]
    names(abundances) <- rownames(my_table)
    abundances_psc_increased <- sum(abundances[psc_increased])
    abundances_psc_decreased <- sum(abundances[psc_decreased])
    
    dys_index <- abundances_psc_increased - abundances_psc_decreased
    dysbiosis_data <- rbind(dysbiosis_data,data.frame(SampleID,PatientID,dys_index))
  }
  colnames(dysbiosis_data) <- c("SampleID","PatientID",name)
  return(dysbiosis_data)
}

## Filtering functions ----

filtering_steps <- function(asv_tab,taxa_tab,metadata,
                            seq_depth_threshold=10000){
  # filtering
  ## sequencing depth 10 000
  data_filt <- seq_depth_filtering(asv_tab,taxa_tab,metadata)
  filt_asv_tab <- data_filt[[1]]
  filt_taxa_tab <- data_filt[[2]]
  filt_metadata <- data_filt[[3]]
  
  seq_step <- dim(filt_asv_tab)[1]
  
  ## NearZeroVar
  data_filt <- nearzerovar_filtering(filt_asv_tab, filt_taxa_tab, filt_metadata)
  filt_asv_tab <- data_filt[[1]]
  filt_taxa_tab <- data_filt[[2]]
  
  nearzero_step <- dim(filt_asv_tab)[1]
  
  # Zerodepth
  # filter the samples with remaining counts smaller than 100
  #colnames(filt_asv_tab[,-1])[colSums(filt_asv_tab[,-1])<100]
  if (TRUE %in% (colSums(filt_asv_tab[,-1])<100)){
    where <- grep(colnames(filt_asv_tab[,-1])[colSums(filt_asv_tab[,-1])<100],colnames(filt_asv_tab))
    filt_asv_tab <- filt_asv_tab[,-where]
    filt_metadata <- filt_metadata[filt_metadata$SampleID %in% colnames(filt_asv_tab),]
    rownames(filt_metadata) <- NULL
    data_checked <- data_check(filt_asv_tab,filt_taxa_tab)
    filt_asv_tab <- data_checked[[1]]
    filt_taxa_tab <- data_checked[[2]]
  }
  
  return(list(filt_asv_tab,filt_taxa_tab, filt_metadata))
}

seq_depth_filtering <- function(asv_tab, taxa_tab, metadata, seq_depth_threshold=10000){
  counts <- apply(asv_tab[,-1],2,sum)
  
  filt_asv_tab <- asv_tab[,c(TRUE,counts>seq_depth_threshold)]
  filt_metadata <- metadata[metadata$SampleID %in% colnames(filt_asv_tab),]
  
  rownames(filt_asv_tab) <- NULL
  rownames(filt_metadata) <- NULL
  
  data_checked <- data_check(filt_asv_tab,taxa_tab)
  filt_asv_tab <- data_checked[[1]]
  filt_taxa_tab <- data_checked[[2]]
  
  return(list(filt_asv_tab, filt_taxa_tab, filt_metadata))
}

prevalence_filtering <- function(prevalence_threshold=0.05,asv_table, taxa_tab, metadata){
  
  groups <- unique(metadata$Group)
  
  filtered_asv_table <- NULL
  asv_table %<>% column_to_rownames("SeqID")
  asvs_to_keep <- c()
  for (group in groups){
    sub_asv_table <- asv_table[,metadata[metadata$Group==group,"SampleID"]]
    sub_asv_table_binary <- sub_asv_table
    sub_asv_table_binary[sub_asv_table>0] <- 1
    prevalence_df <- apply(sub_asv_table_binary,1,sum)/sum(metadata$Group==group)
    prevalence_df <- data.frame(SeqID=rownames(asv_table),prevalence_df)
    
    sub_asv_table <- data.frame(SeqID=rownames(asv_table),sub_asv_table,check.names = FALSE)
    sub_asv_table <- sub_asv_table[prevalence_df$prevalence_df>=prevalence_threshold,]
    asvs_to_keep <- c(asvs_to_keep,rownames(sub_asv_table))
  }
  filtered_asv_table <- asv_table[unique(asvs_to_keep),] %>% rownames_to_column("SeqID")
  data_checked <- data_check(filtered_asv_table,taxa_tab)
  filt_asv_table <- data_checked[[1]]
  filt_taxa_tab <- data_checked[[2]]
  
  
  return(list(filt_asv_table,filt_taxa_tab))
}

abundance_filtering <- function(abundance_threshold=0.05,asv_table, taxa_tab){
  asv_table %<>% column_to_rownames("SeqID")
  
  for (sample in colnames(asv_table)){
    relative_abundances <- asv_table[,sample]/sum(asv_table[,sample])
    where_0 <- relative_abundances < abundance_threshold
    asv_table[where_0,sample] <- 0
  }
  asv_table %<>% rownames_to_column("SeqID")
  data_checked <- data_check(asv_table,taxa_tab)
  filt_asv_table <- data_checked[[1]]
  filt_taxa_tab <- data_checked[[2]]
  
  
  return(list(filt_asv_table,filt_taxa_tab))
}

nearzerovar_filtering <- function(asv_tab,taxa_tab,metadata){
  asv_tab %<>% column_to_rownames("SeqID")
  groups <- unique(metadata$Group)
  
  filtered_asv_table <- NULL
  
  asvs_to_keep <- c()
  for (group in groups){
    sub_asv_table <- asv_tab[,metadata[metadata$Group==group,"SampleID"]]
    
    tab <- sub_asv_table %>% t()
    near_zero <- nearZeroVar(tab)
    asvs_to_keep <- c(asvs_to_keep,rownames(asv_tab[-near_zero,]))
  }
  filtered_asv_table <- asv_tab[unique(asvs_to_keep),] %>% rownames_to_column("SeqID")
  data_checked <- data_check(filtered_asv_table,taxa_tab)
  filt_asv_table <- data_checked[[1]]
  filt_taxa_tab <- data_checked[[2]]
  
  return(list(filt_asv_table,filt_taxa_tab))
}

check_zero_depth <- function(asv_table, taxa_table, metadata){
  if (TRUE %in% (colSums(asv_table[,-1])<100)){
    where <- grep(colnames(asv_table[,-1])[colSums(asv_table[,-1])<100],colnames(asv_table))
    filt_asv_table <- asv_table[,-where]
    filt_metadata <- metadata[metadata$SampleID %in% colnames(filt_asv_table),]
    rownames(filt_metadata) <- NULL
    data_checked <- data_check(filt_asv_table,taxa_table)
    filt_asv_table <- data_checked[[1]]
    filt_taxa_tab <- data_checked[[2]]
    return(list(filt_asv_table, filt_taxa_tab, filt_metadata))
  }
  else (return(list(asv_table,taxa_table,metadata)))
}

final_counts_filtering <- function(asv_tab,filt_asv_tab,filt_metadata,
                                   seq_step,prev_step,nearzero_step){
  df <- data.frame(`Raw data: ASVs`=dim(asv_tab)[1],
                   `Raw data: Samples`=dim(asv_tab)[2]-1,
                   `Sequencing depth filt: ASVs`=seq_step,
                   `Prevalence filt: ASVs` = prev_step,
                   `NearZeroVar filt: ASVs`=nearzero_step,
                   `Filt data: ASVs`=dim(filt_asv_tab)[1],
                   `Filt data: Samples`=dim(filt_asv_tab)[2]-1,
                   `Filt data: Patients`=length(unique(filt_metadata$Patient)),
                   `Filt data: Patients`=length(unique(filt_metadata$subjectid)),
                   `Filtered samples`=ncol(asv_tab[,-1])-ncol(filt_asv_tab[,-1]),
                   Matrices=paste(unique(filt_metadata$Matrix), collapse = ";"),
                   check.names = FALSE)
  
  groups <- c("healthy","non-rPSC","rPSC","pre_ltx","post_ltx","ETOH")
  for (group in groups){
    df[,group] <- sum(filt_metadata$Group==group)
  }
  
  df %<>% t() %>% as.data.frame() 
  
  return(df)
}


## Statistical testing ----

pairwise.wilcox <- function(formula,factors,data, p.adjust.m ='BH')
{
  set.seed(123)
  #co <- combn(unique(as.character(factors)),2)
  co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","non-rPSC","rPSC"),
                   "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","pre_ltx","pre_ltx")) %>% t()
  p_values <- c()
  groups <- c()
  for(elem in 1:ncol(co)){
    x_sub <- data[factors %in% c(co[1,elem],co[2,elem]),]
    group <- paste(co[1,elem],"vs",co[2,elem])
    p_value <- wilcox.test(as.formula(formula), data = x_sub, exact = FALSE)$p.value
    groups <- c(groups,group)
    p_values <- c(p_values,p_value)
  }
  p.adjusted <- p.adjust(p_values,method=p.adjust.m)
  
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  
  df <- data.frame(pair=groups,
                   p.value=p_values,
                   padjust=p.adjusted,
                   sig=sig)
  return(df)

}

lm.model <- function(formula,data){
  model <- coef(summary(lm(as.formula(formula), data = data)))
  p_values <- model[,4]
  sig = c(rep('',length(p_values)))
  sig[p_values <= 0.05] <-'*'
  sig[p_values <= 0.01] <-'**'
  sig[p_values <= 0.001] <-'***'
  df <- data.frame(model,sig=sig)
  return(df)
}

lmer.model <- function(formula,data){
  model <- coef(summary(lmer(as.formula(formula), data = data)))
  p_values <- model[,4]
  sig = c(rep('',length(p_values)))
  sig[p_values <= 0.05] <-'*'
  sig[p_values <= 0.01] <-'**'
  sig[p_values <= 0.001] <-'***'
  df <- data.frame(model,sig=sig)
  return(df)
}


pairwise.lm <- function(formula,factors,data, p.adjust.m ='BH')
{
  set.seed(123)
  #co <- combn(unique(as.character(factors)),2)
  co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","non-rPSC","rPSC"),
                   "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","pre_ltx","pre_ltx")) %>% t()
  
  if ("post_ltx" %in% factors) {
    co <- data.frame("1"=c("healthy","pre_ltx","healthy"),
                     "2"=c("pre_ltx","post_ltx","post_ltx")) %>% t()
  } 
  if (!("pre_ltx" %in% factors)) {
    co <- data.frame("1"=c("non-rPSC","healthy","healthy"),
                     "2"=c("rPSC","rPSC","non-rPSC")) %>% t()
  } 
  if (!FALSE %in% (c("healthy","non-rPSC","rPSC","pre_ltx") %in% factors)){
    co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","pre_ltx","pre_ltx"),
                     "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","non-rPSC","rPSC")) %>% t()
  }
  if ("ibd" %in% factors) {
    co <- data.frame("1"=c("no_ibd"),
                     "2"=c("ibd")) %>% t()
  } 
  models <- c()
  names <- c()
  emeans_models <- c()
  means <- c()
  for(elem in 1:ncol(co)){
    x_sub <- data[factors %in% c(co[1,elem],co[2,elem]),]
    x_sub$Group <- factor(x_sub$Group)
    x_sub$Group <- relevel(x_sub$Group,co[1,elem])
    group <- paste(co[1,elem],"vs",co[2,elem])
    model <- coef(summary(lm(as.formula(formula), data = x_sub)))
    model <- glm_renaming(model,c(co[1,elem],co[2,elem]))
    mean_group <- model[1,]
    if (nrow(model)>2){
      if (model[4,"Pr(>|t|)"]<0.1){
        model_raw <- lm(as.formula(formula), data = x_sub)
        emeans_model <- as.data.frame(contrast(emmeans(model_raw, pairwise ~ Group | Country),method="revpairwise"))
        emeans_models <- rbind(emeans_models,emeans_model)
      } 
    }
    names <- c(names,rownames(model))
    models <- rbind(models,model)
    means <- rbind(mean_group,means)
  }
  
  rownames(means) <- co[1,]
  models <- models[-grep("^[(].+[)]$", names),]
  models %<>% as.data.frame() %>% `row.names<-`(names[-grep("^[(].+[)]$", names)])
  
  if (ncol(co)>1){
    p.adjusted <- p.adjust(models[,"Pr(>|t|)"],method=p.adjust.m)
    models$p.adj <- p.adjusted
  } else models$p.adj <- p.adjusted <- models[,"Pr(>|t|)"]
  
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  
  df <- data.frame(models,
                   sig=sig)
  return(list(df,emeans_models, means))
  
}

pairwise.lmer <- function(formula,factors,data, p.adjust.m ='BH')
{
  set.seed(123)
  #co <- combn(unique(as.character(factors)),2)
  co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","non-rPSC","rPSC"),
                   "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","pre_ltx","pre_ltx")) %>% t()
  
  if ("post_ltx" %in% factors) {
    co <- data.frame("1"=c("healthy","pre_ltx","healthy"),
                     "2"=c("pre_ltx","post_ltx","post_ltx")) %>% t()
  } 
  if (!("pre_ltx" %in% factors)) {
    co <- data.frame("1"=c("non-rPSC","healthy","healthy"),
                     "2"=c("rPSC","rPSC","non-rPSC")) %>% t()
  } 
  if (!FALSE %in% (c("healthy","non-rPSC","rPSC","pre_ltx") %in% factors)){
    co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","pre_ltx","pre_ltx"),
                     "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","non-rPSC","rPSC")) %>% t()
  } 
  if ("ibd" %in% factors) {
    co <- data.frame("1"=c("no_ibd"),
                     "2"=c("ibd")) %>% t()
  } 
  
  models <- c()
  names <- c()
  result_lists <- c()
  means <- c()
  for(elem in 1:ncol(co)){
    x_sub <- data[factors %in% c(co[1,elem],co[2,elem]),]
    x_sub$Group <- factor(x_sub$Group)
    x_sub$Group <- relevel(x_sub$Group,co[1,elem])
    group <- paste(co[1,elem],"vs",co[2,elem])
    model <- coef(summary(lmer(as.formula(formula), data = x_sub)))
    model <- lmer_renaming(model,c(co[1,elem],co[2,elem]))
    mean_group <- model[1,]
    if (nrow(model)>2){
      if (model[4,"Pr(>|t|)"]<0.15){
      # CZ
      group1_group2_cz <- as.data.frame(coef(summary(lmer(as.formula(gsub(" \\* Country","",formula)), data = subset(x_sub,Country=="CZ")))))
      
      # NO
      group1_group2_no <- as.data.frame(coef(summary(lmer(as.formula(gsub(" \\* Country","",formula)), data = subset(x_sub,Country=="NO")))))
      
      # GROUP 1 
      group1_country <- as.data.frame(coef(summary(lmer(as.formula(gsub(" Group \\*","",formula)), data = subset(x_sub,Group==co[1,elem])))))
      
      # GROUP 2 
      group2_country <- as.data.frame(coef(summary(lmer(as.formula(gsub(" Group \\*","",formula)), data = subset(x_sub,Group==co[2,elem])))))
      
      p_values <- c(group1_group2_cz$`Pr(>|t|)`[2], 
                    group1_group2_no$`Pr(>|t|)`[2], 
                    group1_country$`Pr(>|t|)`[2], 
                    group2_country$`Pr(>|t|)`[2])
      
      # Apply Benjamini-Hochberg correction
      p_values_adj <- p.adjust(p_values, method = "BH")
      
      group1_group2_cz$`p.adjusted` <- NA
      group1_group2_no$`p.adjusted` <- NA
      group1_country$`p.adjusted` <- NA
      group2_country$`p.adjusted` <- NA
      
      group1_group2_cz$`p.adjusted`[2] <- p_values_adj[1]
      group1_group2_no$`p.adjusted`[2] <- p_values_adj[2]
      group1_country$`p.adjusted`[2] <- p_values_adj[3]
      group2_country$`p.adjusted`[2] <- p_values_adj[4]
      
      result_list <- list(group1_group2_cz,group1_group2_no,group1_country,group2_country)
      names(result_list) <- c(
        paste(co[1,elem],co[2,elem],"CZ",sep = "_"),
        paste(co[1,elem],co[2,elem],"NO",sep = "_"),
        paste(co[1,elem],"CZ_vs_NO",sep = "_"),
        paste(co[2,elem],"CZ_vs_NO",sep = "_")
      )
      result_lists[[length(result_lists)+1]] <- result_list
      }
    }
    names <- c(names,rownames(model))
    models <- rbind(models,model)
    means <- rbind(mean_group,means)
  }
  rownames(means) <- co[1,]
  models <- models[-grep("^[(].+[)]$", names),]
  models %<>% as.data.frame() %>% `row.names<-`(names[-grep("^[(].+[)]$", names)])
  
  if (ncol(co)>1){
    p.adjusted <- p.adjust(models[,"Pr(>|t|)"],method=p.adjust.m)
    models$p.adj <- p.adjusted
  } else models$p.adj <- p.adjusted <- models[,"Pr(>|t|)"]
  

  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  
  df <- data.frame(models,
                   sig=sig)
  return(list(df,result_lists,means))
  
}

run.adonis <- function(x,factors,covariate=NULL,interaction=FALSE, patients=NULL,
                       sim.function = 'vegdist', sim.method = 'robust.aitchison', 
                       p.adjust.m ='BH',reduce=NULL,perm=999){
  
  x1 = vegdist(x,method=sim.method)
  x2 = data.frame(Fac = factors)
  
  if (!(is.null(covariate))) {
    x2$covariate <- covariate
  }
  
  if (!(is.null(patients))){
    x2$Patient <- patients
    perm <- custom_permutations(x,factors = x2$Fac,patients = x2$Patient)
  }
  
  if (is.null(covariate)){
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm, by="margin");
  } else if ((!is.null(covariate)) & (interaction==FALSE)) {
    ad <- adonis2(x1 ~ Fac + covariate, data = x2,
                  permutations = perm, by="margin");
  } else {
    ad <- adonis2(x1 ~ Fac * covariate, data = x2,
                  permutations = perm, by="terms");
  }
  
  rownames(ad)[1] <- "Group"
  if ("covariate" %in% rownames(ad)) rownames(ad)[2] <- "Country"
  if ("Fac:covariate" %in% rownames(ad)) rownames(ad)[3] <- "Interaction"

  return(ad)
}



pairwise.adonis <- function(x,factors, covariate=NULL, interaction=FALSE,patients=NULL, 
                            sim.function = 'vegdist', sim.method = 'robust.aitchison',
                            p.adjust.m ='BH',reduce=NULL,perm=999)
{
  set.seed(123)
  #co <- combn(unique(as.character(factors)),2)
  co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","non-rPSC","rPSC"),
                   "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","pre_ltx","pre_ltx")) %>% t()
  
  if ("ETOH" %in% factors){
    co <- data.frame("1"=c("healthy","healthy"),
                     "2"=c("ETOH","post_ltx")) %>% t()
  }
  if ("post_ltx" %in% factors) {
    co <- data.frame("1"=c("pre_ltx","pre_ltx","post_ltx"),
                     "2"=c("healthy","post_ltx","healthy")) %>% t()
  } 
  
  if (!("pre_ltx" %in% factors)) {
    co <- data.frame("1"=c("rPSC","rPSC","non-rPSC"),
                     "2"=c("non-rPSC","healthy","healthy")) %>% t()
  } 
  
  if ("ibd" %in% factors) {
    co <- data.frame("1"=c("no_ibd"),
                     "2"=c("ibd")) %>% t()
  } 

  pairs_factor <- c()
  Df_factor <- c()
  SumsOfSqs_factor <- c()
  F.Model_factor <- c()
  R2_factor <- c()
  p.value_factor <- c()
  
  pairs_cov <- c()
  Df_cov <- c()
  SumsOfSqs_cov <- c()
  F.Model_cov <- c()
  R2_cov <- c()
  p.value_cov <- c()
  
  pairs_fac.cov <- c()
  Df_fac.cov <- c()
  SumsOfSqs_fac.cov <- c()
  F.Model_fac.cov <- c()
  R2_fac.cov <- c()
  p.value_fac.cov <- c()
  
  for(elem in 1:ncol(co)){
    x_sub <- x[factors %in% c(co[1,elem],co[2,elem]),]
    x1 = vegdist(x_sub,method=sim.method)
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    if (!(is.null(covariate))) {
      x2$covariate <- covariate[factors %in% c(co[1,elem],co[2,elem])]
    }
    
    if (!(is.null(patients))){
      x2$Patient <- patients[factors %in% c(co[1,elem],co[2,elem])]
      perm <- custom_permutations(x_sub,factors = x2$Fac,patients = x2$Patient)
    }
    
    if (is.null(covariate)){
      ad <- adonis2(x1 ~ Fac, data = x2,
                    permutations = perm, by="margin");
    } else if ((!is.null(covariate)) & (interaction==FALSE)) {
      ad <- adonis2(x1 ~ Fac + covariate, data = x2,
                    permutations = perm, by="margin");
    } else {
      ad <- adonis2(x1 ~ Fac * covariate, data = x2,
                    permutations = perm, by="terms");
    }
    
    pairs_factor <- c(pairs_factor,paste(co[1,elem],'vs',co[2,elem]));
    Df_factor <- c(Df_factor,ad$Df[1])
    SumsOfSqs_factor <- c(SumsOfSqs_factor,ad$SumOfSqs[1])
    F.Model_factor <- c(F.Model_factor,ad$F[1]);
    R2_factor <- c(R2_factor,ad$R2[1]);
    p.value_factor <- c(p.value_factor,ad$`Pr(>F)`[1])
    
    if (!(is.null(covariate))){
      pairs_cov <- c(pairs_cov,paste(co[1,elem],'vs',co[2,elem],", Country"));
      Df_cov <- c(Df_cov,ad$Df[2])
      SumsOfSqs_cov <- c(SumsOfSqs_cov,ad$SumOfSqs[2])
      F.Model_cov <- c(F.Model_cov,ad$F[2]);
      R2_cov <- c(R2_cov,ad$R2[2]);
      p.value_cov <- c(p.value_cov,ad$`Pr(>F)`[2])
      
      if (interaction){
        pairs_fac.cov <- c(pairs_fac.cov,paste(co[1,elem],'vs',co[2,elem],": Country"));
        Df_fac.cov <- c(Df_fac.cov,ad$Df[3])
        SumsOfSqs_fac.cov <- c(SumsOfSqs_fac.cov,ad$SumOfSqs[3])
        F.Model_fac.cov <- c(F.Model_fac.cov,ad$F[3]);
        R2_fac.cov <- c(R2_fac.cov,ad$R2[3]);
        p.value_fac.cov <- c(p.value_fac.cov,ad$`Pr(>F)`[3])
      }
      
    }
    
  }
  p.adjusted_factor <- p.adjust(p.value_factor,method=p.adjust.m)
  if (!(is.null(covariate))) p.adjusted_cov <- p.adjust(p.value_cov,method=p.adjust.m)
  if (!(is.null(covariate))) p.adjusted_fac.cov <- p.adjust(p.value_fac.cov,method=p.adjust.m)
  
  sig_factor = c(rep('',length(p.adjusted_factor)))
  sig_factor[p.adjusted_factor <= 0.05] <-'*'
  sig_factor[p.adjusted_factor <= 0.01] <-'**'
  sig_factor[p.adjusted_factor <= 0.001] <-'***'
  
  if (!(is.null(covariate))) {
    sig_cov = c(rep('',length(p.adjusted_cov)))
    sig_cov[p.adjusted_cov <= 0.05] <-'*'
    sig_cov[p.adjusted_cov <= 0.01] <-'**'
    sig_cov[p.adjusted_cov <= 0.001] <-'***'
  }
  
  if (!(is.null(covariate))) {
    sig_fac.cov = c(rep('',length(p.adjusted_fac.cov)))
    sig_fac.cov[p.adjusted_fac.cov <= 0.05] <-'*'
    sig_fac.cov[p.adjusted_fac.cov <= 0.01] <-'**'
    sig_fac.cov[p.adjusted_fac.cov <= 0.001] <-'***'
  }
  
  pairw.res_factor <- data.frame(pairs_factor,Df_factor,SumsOfSqs_factor,F.Model_factor,R2_factor,p.value_factor,p.adjusted_factor,sig_factor)
  if (!(is.null(covariate))) pairw.res_cov <- data.frame(pairs_cov,Df_cov,SumsOfSqs_cov,F.Model_cov,R2_cov,p.value_cov,p.adjusted_cov,sig_cov)
  if (!(is.null(covariate))) pairw.res_fac.cov <- data.frame(pairs_fac.cov,Df_fac.cov,SumsOfSqs_fac.cov,F.Model_fac.cov,R2_fac.cov,p.value_fac.cov,p.adjusted_fac.cov,sig_fac.cov)
  
  class(pairw.res_factor) <- c("pwadonis", "data.frame")
  if (!(is.null(covariate))) class(pairw.res_cov) <- c("pwadonis", "data.frame")
  if (!(is.null(covariate))) class(pairw.res_fac.cov) <- c("pwadonis", "data.frame")
  
  if (!(is.null(covariate))) return(list(pairw.res_factor,pairw.res_cov,pairw.res_fac.cov))
  else (return(pairw.res_factor))
}

adonis_postanalysis <- function(x,factors, covariate, group1, group2,
                                patients = NULL, 
                                sim.method = 'robust.aitchison'){
  
  # CZ
  x_sub <- x[(factors %in% c(group1,group2)) & covariate=="CZ",]
  
  x1 = vegdist(x_sub,method=sim.method)
  x2 = data.frame(Fac = factors[(factors %in% c(group1,group2)) & covariate=="CZ"])
  
  if (!(is.null(patients))){
    x2$Patient <- patients[(factors %in% c(group1,group2)) & covariate=="CZ"]
    perm <- custom_permutations(x_sub,factors = x2$Fac,patients = x2$Patient)
  } else {
    perm <- 999
  }
  
  group1_group2_cz <- adonis2(x1 ~ Fac, data = x2,
                              permutations = perm)
  
  # NO
  x_sub <- x[(factors %in% c(group1,group2)) & covariate=="NO",]
  
  x1 = vegdist(x_sub,method=sim.method)
  x2 = data.frame(Fac = factors[(factors %in% c(group1,group2)) & covariate=="NO"])
  
  if (!(is.null(patients))){
    x2$Patient <- patients[(factors %in% c(group1,group2)) & covariate=="NO"]
    perm <- custom_permutations(x_sub,factors = x2$Fac,patients = x2$Patient)
  } else {
    perm <- 999
  }
  
  group1_group2_no <- adonis2(x1 ~ Fac, data = x2,
                              permutations = perm)
  
  # group1
  x_sub <- x[(factors %in% c(group1)),]
  x1 = vegdist(x_sub,method=sim.method)
  x2 = data.frame(Fac = factors[(factors %in% c(group1))],
                  Cov = covariate[(factors %in% c(group1))])
  
  if (!(is.null(patients))){
    x2$Patient <- patients[(factors %in% c(group1))]
    perm <- custom_permutations(x_sub,factors = x2$Fac,patients = x2$Patient)
  } else {
    perm <- 999
  }
  
  group1_country <- adonis2(x1 ~ Cov, data = x2,
                            permutations = perm)
  
  # group2
  x_sub <- x[(factors %in% c(group2)),]
  x1 = vegdist(x_sub,method=sim.method)
  x2 = data.frame(Fac = factors[(factors %in% c(group2))],
                  Cov = covariate[(factors %in% c(group2))])
  
  if (!(is.null(patients))){
    x2$Patient <- patients[(factors %in% c(group2))]
    perm <- custom_permutations(x_sub,factors = x2$Fac,patients = x2$Patient)
  } else {
    perm <- 999
  }
  
  
  group2_country <- adonis2(x1 ~ Cov, data = x2,
                            permutations = perm)
  
  p_values <- c(group1_group2_cz$`Pr(>F)`[1], 
                group1_group2_no$`Pr(>F)`[1], 
                group1_country$`Pr(>F)`[1], group2_country$`Pr(>F)`[1])
  
  # Apply Benjamini-Hochberg correction
  p_values_adj <- p.adjust(p_values, method = "BH")
  
  group1_group2_cz$`p.adjusted` <- NA
  group1_group2_no$`p.adjusted` <- NA
  group1_country$`p.adjusted` <- NA
  group2_country$`p.adjusted` <- NA
  
  group1_group2_cz$`p.adjusted`[1] <- p_values_adj[1]
  group1_group2_no$`p.adjusted`[1] <- p_values_adj[2]
  group1_country$`p.adjusted`[1] <- p_values_adj[3]
  group2_country$`p.adjusted`[1] <- p_values_adj[4]
  
  result_list <- list(group1_group2_cz,group1_group2_no,group1_country,group2_country)
  names(result_list) <- c(
    paste(group1,group2,"CZ",sep = "_"),
    paste(group1,group2,"NO",sep = "_"),
    paste(group1,"CZ_vs_NO",sep = "_"),
    paste(group2,"CZ_vs_NO",sep = "_")
  )
  return(result_list)
}

pairwise.adonis_A <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=9999)
{
  # This function performs pairwise PERMANOVA, code is modified version of 
  # https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
  # inputs:
  # x - asv table
  # factors - vector with groups to be tested
  # sim.function - function for distance calculation
  # sim.method - distance metric
  # p.adjust.m - FDR method
  # perm - number of permutations
  # outputs:
  # data frame with pairwise PERMANOVA results
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

### Method summary
summary.pwadonis = function(object, ...) {
  cat("Result of pp <- ordering_pp(pairwise.adonis:\n")
  cat("\n")
  print(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}


custom_permutations <- function(x,factors,patients, perm=999){
  # factors - vector of the labels
  # patients - vector of the pacients
  # perm - number of permutations
  
  set.seed(123)
  # create data frame
  patients_df <- data.frame(patient=patients,
                            group=factors)
  # get unique patients
  patients_df <- patients_df[!duplicated(patients_df$patient),]
  
  # permutations across patients
  perm <- how(nperm = perm)
  perm <- shuffleSet(nrow(patients_df), control = perm)
  colnames(perm) <- patients_df$patient
  
  perm_final <- matrix(nrow=nrow(perm),ncol = length(factors)) 
  colnames(perm_final) <- rownames(x)
  
  # one patient = same permutation index
  for (j in 1:ncol(perm_final)){
    which_patient <- patients[j]
    which_ind_patient <- perm[,as.character(which_patient)]
    perm_final[,j] <- which_ind_patient
  }
  return(perm_final)
}

ordering_pp <- function(pp){
  groups_order <- c("healthy vs non-rPSC","healthy vs rPSC","healthy vs pre_ltx","non-rPSC vs rPSC","non-rPSC vs pre_ltx","rPSC vs pre_ltx")
  
  ordering <- c()
  for (group_pair in groups_order) {
    where <- grep(paste0("^",group_pair), pp$pairs)
    
    if (length(where)==0) {
      flip <- unlist(strsplit(group_pair,split=" vs "))
      flip <- paste0(flip[2]," vs ",flip[1])
      where <- grep(paste0("^",flip), pp$pairs)
    }
    ordering <- c(ordering,where)
  }
  
  pp <- pp[ordering,]
  return(pp)
  
}


replace_p_values <- function(p) {
  if (p <= 0.001) {
    return("***")
  } else if (p <= 0.01) {
    return("**")
  } else if (p <= 0.05) {
    return("*")
  } else {
    return(as.character(p))
  }
}

add_significance <- function(model_df){
  model_df$sig <- c(rep('',nrow(model_df)))
  model_df$sig[model_df$`Pr(>|t|)` <= 0.05] <-'**'
  model_df$sig[model_df$`Pr(>|t|)` <= 0.01] <-'**'
  model_df$sig[model_df$`Pr(>|t|)` <= 0.001] <-'***'
  return(model_df)
}

glm_renaming <- function(glm_data, group){
  # This function renames the results of linDA (linda() function $output) based on the chosen group[1] 
  # inputs:
  # linda_data - result of linda()$output
  # group[1] - strings corresponding to the name of the group[1]
  # outputs:
  # linda_data - renamed data
  
  if (is.data.frame(glm_data) | is.matrix(glm_data)){
    # renaming group[1]
    rownames(glm_data) <- gsub("Intercept",group[1], rownames(glm_data))
    
    # renaming main effects
    rownames(glm_data)[grepl("Group",rownames(glm_data)) & !(grepl(":",rownames(glm_data)))] <-  paste(group[1],"vs",rownames(glm_data)[grepl("Group",rownames(glm_data)) & !(grepl(":",rownames(glm_data)))])
    rownames(glm_data)[grepl("Country",rownames(glm_data)) & !(grepl(":",rownames(glm_data)))] <- paste(group[1], ",", group[2], "-", "CZ vs NO")
    
    # renaming interactions
    rownames(glm_data)[grepl(":",rownames(glm_data))] <- paste(group[1],"vs",rownames(glm_data)[grepl(":",rownames(glm_data))])
    
  } else{
    message("Check the glm data class")
  }
  return(glm_data)
}

lmer_renaming <- function(lmer_data, group){
  # This function renames the results of linDA (linda() function $output) based on the chosen group[1] 
  # inputs:
  # linda_data - result of linda()$output
  # group[1] - strings corresponding to the name of the group[1]
  # outputs:
  # linda_data - renamed data
  
  if (is.data.frame(lmer_data) | is.matrix(lmer_data)){
    # renaming group[1]
    rownames(lmer_data) <- gsub("intercept",group[1], rownames(lmer_data))
    
    # renaming main effects
    rownames(lmer_data)[grepl("Group",rownames(lmer_data)) & !(grepl(":",rownames(lmer_data)))] <-  paste(group[1],"vs",rownames(lmer_data)[grepl("Group",rownames(lmer_data)) & !(grepl(":",rownames(lmer_data)))])
    rownames(lmer_data)[grepl("Country",rownames(lmer_data)) & !(grepl(":",rownames(lmer_data)))] <- paste(group[1],"vs",group[2], "- CZ vs NO")
    #rownames(lmer_data)[grepl("Matrix",rownames(lmer_data)) & !(grepl(":",rownames(lmer_data)))] <- paste(group[1]," vs ",group[2], "Cecum vs", rownames(lmer_data)[grepl("Matrix",rownames(lmer_data)) & !(grepl(":",rownames(lmer_data)))])
    
    # renaming interactions
    rownames(lmer_data)[grepl(":",rownames(lmer_data))] <- paste(group[1],"vs",rownames(lmer_data)[grepl(":",rownames(lmer_data))])
    
  } else{
    message("Check the lmer data class")
  }
  return(lmer_data)
}

## Visualization functions ----
horizontal_barplot <- function(wb,taxa){
  sheets_names <- sheets(wb)
  groups_names <- unique(unlist(strsplit(sheets_names,split=" vs ")))
  prevalences_df <- NULL
  
  for (name in sheets_names){
    group1_name <- unlist(strsplit(name,split=" vs "))[1]
    group2_name <- unlist(strsplit(name,split=" vs "))[2]
    df <- readWorkbook(wb, sheet = name)
    #df <- df[df$final_sig,]
    if (is.null(prevalences_df)){
      prevalences_df <- df[,c("SeqID",paste0("PREVALENCE_.",group1_name),
                              paste0("PREVALENCE_.",group2_name))] %>% `rownames<-`(NULL)
    } else{
      prevalences_df <- merge(prevalences_df,
                               df[,c("SeqID",paste0("PREVALENCE_.",group1_name),
                                     paste0("PREVALENCE_.",group2_name))] %>% `rownames<-`(NULL),all=TRUE)
    }
  }
  
  # Clean the dataframe by grouping and filling NAs
  if (any(grepl("post_ltx",colnames(prevalences_df)))){
    prevalences_df <- prevalences_df %>%
    group_by(SeqID) %>%
    summarise(
      PREVALENCE_.healthy = coalesce(PREVALENCE_.healthy[1], PREVALENCE_.healthy[2]),
      PREVALENCE_.post_ltx = coalesce(PREVALENCE_.post_ltx[1], PREVALENCE_.post_ltx[2]),
      PREVALENCE_.pre_ltx = coalesce(PREVALENCE_.pre_ltx[1], PREVALENCE_.pre_ltx[2])
    ) 
  } else if (any(grepl("rPSC",colnames(prevalences_df)))){
    prevalences_df <- prevalences_df %>%
      group_by(SeqID) %>%
      summarise(
        PREVALENCE_.healthy = coalesce(PREVALENCE_.healthy[1], PREVALENCE_.healthy[2]),
        `PREVALENCE_.non-rPSC` = coalesce(`PREVALENCE_.non-rPSC`[1], `PREVALENCE_.non-rPSC`[2]),
        PREVALENCE_.rPSC = coalesce(PREVALENCE_.rPSC[1], PREVALENCE_.rPSC[2])
      ) 
  } else if (any(grepl("ibd",colnames(prevalences_df)))){
    prevalences_df <- prevalences_df %>%
      group_by(SeqID) %>%
      summarise(
        PREVALENCE_.no_ibd = coalesce(PREVALENCE_.no_ibd[1], PREVALENCE_.no_ibd[2]),
        `PREVALENCE_.ibd` = coalesce(`PREVALENCE_.ibd`[1], `PREVALENCE_.ibd`[2]),
      ) 
    } else cat("CHYBA")
  
  

  prevalences_df <- prevalences_df %>% dplyr::filter(SeqID %in% taxa) %>%
    dplyr::mutate(SeqID=factor(SeqID,levels=taxa))
  
  
  prevalences_df_melt <- melt(prevalences_df)
  
  prevalences_df_melt <- prevalences_df_melt %>%
    mutate(variable = case_when(
      variable == "PREVALENCE_.healthy" ~ "HC",
      variable == "PREVALENCE_.post_ltx" ~ "post_LTx",
      variable == "PREVALENCE_.pre_ltx" ~ "pre_LTx",
      variable == "PREVALENCE_.rPSC" ~ "rPSC",
      variable == "PREVALENCE_.non-rPSC" ~ "non-rPSC",
      variable == "PREVALENCE_.ibd" ~ "ibd",
      variable == "PREVALENCE_.no_ibd" ~ "no_ibd",
      TRUE ~ variable  # Keep other values as is
    ))
  
  if ("post_LTx" %in% prevalences_df_melt$variable) {
    colors <- c("#309f87","#f9c675","#425387","#d55c4a")
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("HC","pre_LTx","post_LTx"))
  }
  else if ("rPSC" %in% prevalences_df_melt$variable) {
    colors <- c("#309f87","#F08080","#A00000")
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("HC","non-rPSC","rPSC"))
  } else if ("ibd" %in% alpha_data$Group){
    #colors <- c("#1B7837","#B2182B")  
    #colors <- c("#D04E36", "#3F7D3C")  
    colors <- c("#A06A2C", "#B2182B")  
    
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("no_ibd","ibd"))
  }
  else (cat("chyba","\n"))
  
  p <- ggplot(prevalences_df_melt, aes(x = SeqID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    coord_flip() + 
    scale_fill_manual(values=colors)  +
    theme_minimal() + 
    theme(axis.title.x =element_text(size=10),
          axis.text.x = element_text(size=4),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          axis.line.y.left = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.ticks.y.left = element_blank())+
    ylab("Prevalence") + 
    xlab("") + 
    theme(legend.position = "none") #+ 
    #scale_y_reverse()

  return(p)
}

mpse_barplot <- function(mpse_object, taxonomic_level="Genus", topn=20){
  p <- ikem_mpse %>%
    mp_plot_abundance(
      .abundance = Abundance,
      .group= Group_matrix,
      taxa.class = !!as.name(taxonomic_level), 
      topn = topn,
      relative = TRUE,
      plot.group = TRUE,
      force = TRUE)
  return(p)
}

boxplot_subgroups <-function(alpha_data,labels){
  alpha_data <- alpha_data %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  
  if ("post_LTx" %in% alpha_data$Group) {
    colors <- c("#309f87","#f9c675","#425387","#d55c4a")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","pre_LTx","post_LTx"))
  }
  else if ("rPSC" %in% alpha_data$Group) {
    colors <- c("#309f87","#F08080","#A00000")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","non-rPSC","rPSC"))
  }
  else (cat("chyba","\n"))
  alpha_data_melted <- melt(alpha_data)
  
  
  p <- ggplot() + 
    geom_boxplot(data=alpha_data_melted, aes(x=Group, y=value, fill=Group, color=Country)) + 
    geom_text(data = labels, aes(x=x, y= y, label = label),size=3) +
    scale_color_manual(values=c("#000000","#000000","#000000","#000000")) +
    scale_fill_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) + 
    guides(fill="none",color = "none") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    facet_wrap(~variable, ncol = 4,scales = "free")
  return(p)
}

alpha_diversity_countries <-function(alpha_data,show_legend=FALSE){
  alpha_data <- alpha_data %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  
  if ("post_LTx" %in% alpha_data$Group) {
    colors <- c("#309f87","#f9c675","#425387","#d55c4a")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","pre_LTx","post_LTx"))
  } else if ("rPSC" %in% alpha_data$Group) {
    colors <- c("#309f87","#F08080","#A00000")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","non-rPSC","rPSC"))
  } else if ("ibd" %in% alpha_data$Group){
    #colors <- c("#1B7837","#B2182B")  
    #colors <- c("#D04E36", "#3F7D3C")  
    colors <- c("#A06A2C", "#B2182B")  

    alpha_data$Group <- factor(alpha_data$Group,levels = c("no_ibd","ibd"))
  }
  else (print("chyba","\n"))
  
  alpha_data$Country <- factor(alpha_data$Country,levels = c('CZ','NO'))
  alpha_richness_data <- melt(alpha_data) %>% dplyr::filter(variable=="Observe")
  
  richness_limit <- max(alpha_data$Observe) + 0.2*max(alpha_data$Observe)
  
  p_richness <- ggplot(alpha_richness_data, aes(x=Group, y=value, fill=Country)) + 
    geom_boxplot(aes(fill=Country), outlier.shape = NA, position=position_dodge(width=0.8)) + 
    geom_jitter(aes(x=Group, y=value, color=Group, shape=Country), 
                position=position_jitterdodge(jitter.width=0.7, dodge.width=0.8), 
                size=1) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    xlab("") + 
    ylab("Richness") + 
    scale_y_continuous(breaks = seq(0, richness_limit, by = 50)) + 
    ylim(0,richness_limit)  + 
    theme(axis.text.x = element_text(angle = 45,face = "bold",vjust = 0.5))
  
  alpha_shannon_data <- melt(alpha_data) %>% dplyr::filter(variable=="Shannon")
  shannon_limit <- max(alpha_data$Shannon) + 0.2*max(alpha_data$Shannon)
  
  p_shannon <- ggplot(alpha_shannon_data, aes(x=Group, y=value, fill=Country)) + 
    geom_boxplot(aes(fill=Country), outlier.shape = NA, position=position_dodge(width=0.8)) + 
    geom_jitter(aes(x=Group, y=value, color=Group, shape=Country), 
                position=position_jitterdodge(jitter.width=0.7, dodge.width=0.8), 
                size=1) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    xlab("") + 
    ylab("Shannon") + 
    scale_y_continuous(breaks = seq(0, shannon_limit, by = 50)) + 
    ylim(0,shannon_limit) + 
    theme(axis.text.x = element_text(angle = 45,face = "bold",vjust = 0.5))
  
  if (show_legend){
    p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="right")
  } else {
    p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="none")
  }
  
  return(p)
}

alpha_diversity_custom <- function(alpha_data, size=1.5){
  alpha_data_melted <- melt(alpha_data)
  if ("post_ltx" %in% alpha_data$Group) colors <- c("#309f87","#425387","#f9c675","#d55c4a")
  else if ("rPSC" %in% alpha_data$Group) colors <- c("#309f87","#F08080","#A00000")
  else (cat("chyba","\n"))
  p <- ggplot() + 
    geom_boxplot(data=alpha_data_melted, aes(x=Group, y=value),outliers = FALSE) + 
    geom_jitter(width = 0.2,height = 0,data=alpha_data_melted,aes(x=Group, y=value, color=Group),size=size) +
    #scale_color_manual(values=c("#000000","#000000","#000000","#000000")) +
    scale_color_manual(values=colors) + 
    scale_fill_manual(values=colors) + 
    guides(fill="none",color = "none") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(axis.text.x = element_text(angle = 0,face = "bold"))+ 
    facet_wrap(~variable, ncol = 4,scales = "free")  + 
    xlab("") + 
    ylab("Alpha index value")
  
  return(p)
}

alpha_diversity_custom_2 <- function(alpha_data, size=1.5,width=0.2){
  #colnames(alpha_data)[which(colnames(alpha_data)=="Observe")] <- "Richness"
  alpha_data <- alpha_data %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  
  if ("post_LTx" %in% alpha_data$Group) {
    colors <- c("#309f87","#f9c675","#425387","#d55c4a")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","pre_LTx","post_LTx"))
  }
  else if ("rPSC" %in% alpha_data$Group) {
    colors <- c("#309f87","#F08080","#A00000")
    alpha_data$Group <- factor(alpha_data$Group,levels = c("HC","non-rPSC","rPSC"))
  } else if ("ibd" %in% alpha_data$Group){
    #colors <- c("#1B7837","#B2182B")  
    #colors <- c("#D04E36", "#3F7D3C")  
    colors <- c("#A06A2C", "#B2182B")  
    
    alpha_data$Group <- factor(alpha_data$Group,levels = c("no_ibd","ibd"))
  } else (print("chyba","\n"))
  richness_limit <- max(alpha_data$Richness) + 0.2*max(alpha_data$Richness)
  p_richness <- ggplot() + 
    geom_boxplot(data=alpha_data, aes(x=Group, y=Richness),outliers = FALSE) + 
    geom_jitter(width = width,height = 0,data=alpha_data,aes(x=Group, y=Richness, color=Group,
                                                           shape=Country),size=size) +
    scale_y_continuous(breaks = seq(0, richness_limit, by = 50)) + 
    ylim(0,richness_limit) + 
    #scale_color_manual(values=c("#000000","#000000","#000000","#000000")) +
    scale_color_manual(values=colors) + 
    scale_fill_manual(values=colors) + 
    guides(fill="none",color = "none") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    theme(axis.text.x = element_text(angle = 45,face = "bold",vjust = 0.5))+ 
    xlab("") + 
    ylab("Richness")
  
  shannon_limit <- max(alpha_data$Shannon) + 0.2*max(alpha_data$Shannon)
  p_shannon <- ggplot() + 
    geom_boxplot(data=alpha_data, aes(x=Group, y=Shannon),outliers = FALSE) + 
    geom_jitter(width = width,height = 0,data=alpha_data,aes(x=Group, y=Shannon, color=Group,
                                                           shape=Country),size=size) +
    scale_y_continuous(breaks = seq(0, shannon_limit, by = 1)) + 
    ylim(0,shannon_limit) + 
    #scale_color_manual(values=c("#000000","#000000","#000000","#000000")) +
    scale_color_manual(values=colors) + 
    scale_fill_manual(values=colors) + 
    guides(fill="none",color = "none") + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45,face = "bold",vjust = 0.5))+ 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    xlab("") + 
    ylab("Shannon")
  
  p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="none")
  return(p)
}


volcano_plot <- function(glm.test,glm.eff,group,taxa_table,cutoff.pval=0.05,cutoff.rare=0){
  column <- paste(group,":pval.padj", sep="")
  called <- glm.test[,column] <= cutoff.pval
  all.p <- glm.test[,column]
  all.col <- rgb(0, 0, 0, 0.2)
  
  data_df <- data.frame(x=glm.eff[[group]]$diff.btw,
                        y=-1*log10(all.p),
                        rab.all=glm.eff[[group]]$rab.all,
                        name=rownames(glm.eff[[group]]))
  data_rare <- data_df[data_df$rab.all<cutoff.rare,]
  data_called <- data_df[called,]
  taxa_table <- taxa_table %>% column_to_rownames("SeqID")
  data_called$name <- taxa_table[data_called$name,"Genus"]
  
  maximum <- max(c(abs(min(data_df$x)), abs(max(data_df$x))))
  p <- ggplot(data=data_df, aes(x=x,y=y)) +
    geom_point(colour=all.col, shape=19, size=3) + theme_bw() + 
    xlab("Median Log"[2]~" Difference") + ylab("-1 * Median Log"[10]~" q value") +
    geom_point(data=data_called, aes(x=x,y=y), colour="red", size=3)+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15),
          axis.title=element_text(size=10)) + 
    geom_vline(xintercept=1.5, linetype=2, colour="grey")+ 
    geom_vline(xintercept=-1.5, linetype=2, colour="grey")+ 
    coord_cartesian(xlim = c(-maximum,maximum), clip="on") +
    geom_hline(yintercept=-1*log10(0.05),linetype=2, colour="grey") +
    geom_text_repel(data=data_called, aes(x = x, y = y, 
                                          label = name),size=3)
  return(p)
}

pca_plots <- function(mpse_object,by="Group"){
  plot_12 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    ellipse = TRUE, show.legend = FALSE) +
    scale_fill_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC1 vs PC2")
  
  plot_13 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    .dim=c(1,3),
    ellipse = TRUE, show.legend = FALSE ) +
    scale_fill_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC1 vs PC3")
  
  plot_23 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    .dim=c(2,3),
    ellipse = TRUE, show.legend = FALSE ) +
    scale_fill_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#00A087FF", "#3C5488FF","#FFC470","#DD5746"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC2 vs PC3")
  
  return(list(plot_12,plot_13,plot_23))
  
}

ordiArrowMul_custom <- function (x, ord, at = c(0,0), fill = 0.75,
                            display, choices = c(1,2)) {
  X <- do.call(rbind,scores(x,c("vectors", "factors")))
  u <- c(0,max(abs(range(ord$vectors[,1]))),0,max(abs(range(ord$vectors[,2]))))
  u <- u - rep(at, each = 2)
  r <- c(range(X[,1], na.rm = TRUE), range(X[,2], na.rm = TRUE))
  ## 'rev' takes care of reversed axes like xlim(1,-1)
  rev <- sign(diff(u))[-2]
  if (rev[1] < 0)
    u[1:2] <- u[2:1]
  if (rev[2] < 0)
    u[3:4] <- u[4:3]
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  return(fill * min(u))
}

pca_plot_custom <- function(asv_table,taxa_table,metadata, 
                            measure="robust.aitchison",
                            show_boxplots = TRUE,
                            variable = "Group", size=2,
                            show_legend=TRUE,
                            clinical=FALSE,
                            clinical_metadata=NULL){
  
  metadata <- metadata %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  
  #metadata$Group <- factor(metadata$Group,levels = c("HC","pre_LTx","non-rPSC","rPSC"))
  if (clinical & !is.null(clinical_metadata)){
    clinical_metadata <- clinical_metadata %>%
      mutate(Group = case_when(
        Group == "healthy" ~ "HC",
        Group == "post_ltx" ~ "post_LTx",
        Group == "pre_ltx" ~ "pre_LTx",
        TRUE ~ Group  # Keep other values as is
      )) %>% 
      dplyr::filter(SampleID %in% metadata$SampleID) %>%
      column_to_rownames("SampleID") %>%
      dplyr::select(-c(Matrix,Group,Country))
    
    if (length(unique(clinical_metadata$PatientID))==nrow(clinical_metadata)){
      clinical_metadata %<>% dplyr::select(-c(PatientID))
    }
  }
  
  if ("post_LTx" %in% metadata$Group) {
    colors <- c("#309f87","#f9c675","#425387")
    metadata$Group <- factor(metadata$Group,levels = c("HC","pre_LTx","post_LTx"))
  }
  else if (("rPSC" %in% metadata$Group) & 
           ("pre_LTx" %in% metadata$Group)) {
    colors <- c("#309f87","#f9c675","#F08080","#A00000")
    metadata$Group <- factor(metadata$Group,levels = c("HC","pre_LTx","non-rPSC","rPSC"))
  }
  else if ("ibd" %in% metadata$Group){
    #colors <- c("#1B7837","#B2182B")  
    #colors <- c("#D04E36", "#3F7D3C")  
    colors <- c("#A06A2C", "#B2182B")  
    
    metadata$Group <- factor(metadata$Group,levels = c("no_ibd","ibd"))
  }
  else {colors <- c("#309f87","#F08080","#A00000")
  metadata$Group <- factor(metadata$Group,levels = c("HC","non-rPSC","rPSC"))
  }
  
  ps <- construct_phyloseq(asv_table,taxa_table,metadata)
  if (measure=="robust.aitchison"){
    ps <- microbiome::transform(ps,"rclr")
    measure_final <- "euclidean"
  } else if (measure=="bray") {
    ps <- transform_sample_counts(ps, function(x) x/sum(x))
    measure_final <- "bray"
  } else if (measure=="jaccard") {
    measure_final <- "jaccard"
  }
  
  data_for_pca <- as.data.frame(t(asv_table[,-which(colnames(asv_table)=="SeqID")]))
  data_for_pca <- data_for_pca[metadata$SampleID,]
  data_for_pca <- cbind(data_for_pca,metadata)
  
  distMat <- phyloseq::distance(ps, method = measure_final, type = "samples")
  
  ord <- phyloseq::ordinate(ps, method = "PCoA", distance = distMat)
  
  imp_vec <- ord$values$Relative_eig
  pca_vec <- ord$vectors
  
  x_lab = paste("PCo1 ", "(",round(imp_vec[1]*100,2),"%", ")", sep="")
  y_lab = paste("PCo2 ", "(",round(imp_vec[2]*100,2),"%", ")", sep="")
  
  p <- ggplot(data_for_pca) + 
    geom_point(aes(x=pca_vec[,1],
                   y=pca_vec[,2],
                   color=!!sym(variable),
                   shape=Country),
               show.legend = show_legend,size=size) + 
    stat_ellipse(aes(x=pca_vec[,1],y=pca_vec[,2],color=!!sym(variable)),show.legend = FALSE) + 
    xlab(x_lab)+
    ylab(y_lab)+
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(legend.position = "right") + 
    scale_color_manual(values=colors)
  
  if (clinical){
    set.seed(123)
    if ("PatientID" %in% colnames(clinical_metadata)){
      my_df <- ord$vectors
      clinical_metadata <- clinical_metadata[rownames(my_df),]
      #my_df <- my_df[rownames(clinical_metadata),]
      my_df <- cbind(my_df,PatientID=clinical_metadata$PatientID) %>% 
        as.data.frame() %>% 
        rownames_to_column("SampleID")
      my_df <- my_df %>%
        #group_by(PatientID) %>%
        distinct(PatientID, .keep_all = TRUE) %>%
        as.data.frame() %>%
        column_to_rownames("SampleID") %>%
        dplyr::select(-PatientID) #%>%
        #as.numeric() %>%
        #as.matrix() 
      sample_names <- rownames(my_df)
      ord$vectors <- as.data.frame((lapply(my_df, as.numeric))) %>%
        `rownames<-`(sample_names) %>% as.matrix()
      clinical_metadata %<>% dplyr::select(-PatientID)
      clinical_metadata <- clinical_metadata[rownames(ord$vectors),]
    }
    
    fit <- envfit(ord=ord$vectors, env=clinical_metadata, perm = 999,na.rm=TRUE)
    vector_coordinates <- as.data.frame(scores(fit, "vectors")) * ordiArrowMul_custom(fit, ord)
    vector_coordinates <- vector_coordinates[fit$vectors$pvals < 0.1,]
    factor_coordinates <- as.data.frame(scores(fit, "factors")) * ordiArrowMul_custom(fit,ord)
    factor_coordinates <- factor_coordinates[fit$factors$pvals < 0.1,]
    pca_loadings <- rbind(vector_coordinates,factor_coordinates)
    
    # Radial shift function
    rshift = function(r, theta, a=0.03, b=0.03) {
      r + a + b*abs(cos(theta))
    }
    
    pca_loadings %<>% dplyr::mutate(r = sqrt(Axis.1^2 + Axis.2^2),
                                   theta = atan2(Axis.2,Axis.1),
                                   rnew = rshift(r, theta),
                                   xnew = rnew*cos(theta),
                                   ynew = rnew*sin(theta)) %>%
      rownames_to_column("Variable")
    
    p <- p + 
      geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
                          arrow = arrow(length = unit(0.1, "cm")), 
                          color = "black") + 
      #geom_text(data = pca_loadings, aes(x = xnew, y = ynew, label = Variable), size=2)
    ggrepel::geom_text_repel(data = pca_loadings, aes(x = xnew, y = ynew, label = Variable), 
                             size=2.5, force = 0.01,force_pull = 1)
  }
  
  if (show_boxplots){
    pmain <- p
    xdens <- axis_canvas(pmain,axis="x") +  
      geom_boxplot(data=data_for_pca, aes(x=pca_vec[,1],fill=!!sym(variable))) + 
      scale_fill_manual(values=colors) + 
      expand_limits(y = -0.4)
    ydens = axis_canvas(pmain,axis="y", coord_flip = TRUE) + 
      geom_boxplot(data=data_for_pca,aes(x=pca_vec[,2],fill=!!sym(variable))) + 
      scale_fill_manual(values=colors) + 
      coord_flip() + expand_limits(y = -0.4)
    
    p <- insert_xaxis_grob(pmain,xdens,grid::unit(.2,"null"), position="top",clip="off")
    p <- insert_yaxis_grob(p,ydens,grid::unit(.2,"null"), position="right",clip="off")
    
    p <- as.ggplot(ggdraw(p))
  }
  return(p)
}


umap_plot <- function(asv_table,metadata, distance="robust.aitchison",neighbors=10,
                      min_dist=0.1,by=NULL){
  custom.config <- umap.defaults
  custom.config$random_state <- 123
  custom.config$n_neighbors <- neighbors
  custom.config$min_dist <- min_dist
  
  r.dist <- as.matrix(vegdist(t(asv_table %>% column_to_rownames("SeqID")), method=distance))
  umap_plot <- umap(r.dist,input="dist",config=custom.config)
  data_umap <- data.frame(umap_plot$layout)
  data_umap <- merge(data_umap %>% rownames_to_column("SampleID"),metadata,by="SampleID")
  if ("Country" %in% colnames(data_umap)){
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=Group, shape=Country)) + geom_point() + theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") + 
      scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a"))
      
  } else {
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=Group)) + geom_point() + theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") + 
      scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a"))
  }
  
  return(p)
}

is_dna_sequence <- function(sequence) {
  # This function serves as a control for DNA sequence recognition. 
  # Check if the sequence contains only valid DNA characters (A, T, C, G)
  if (grepl("^[ATCGatcg]+$", sequence)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

heatmap_correlation <- function(corrs){
  
}


heatmap_linda <- function(linda.output,taxa_tab){

  # Create a data frame
  linda.output <- lapply(linda.output, function(df) df[, c("padj", "log2FoldChange")])
  groups <- names(linda.output)
  groups <- rep(groups, each = 2)
  
  # Convert rownames to a column
  linda.output <- map(linda.output, ~ .x %>% rownames_to_column("SeqID"))
  
  # Perform full outer join on the list of dataframes
  linda.output <- purrr::reduce(linda.output, full_join, by = "SeqID")
  
  colnames(linda.output) <- c("SeqID",paste0(groups,c(":padj",":log2FoldChange")))
  
  #names_cols <- linda.output %<>% dplyr::bind_cols(.name_repair = "unique")
  
  or_logic <- rep(FALSE,dim(linda.output)[1])
  for (col in grep(":padj",colnames(linda.output))){
    new_or_logic <- (linda.output[,col] < 0.1)
    new_or_logic[is.na(new_or_logic)] <- FALSE
    or_logic <- (or_logic) | (new_or_logic)
  }
  
  linda.output_reduced <- linda.output[or_logic,]
  
  # Extract relevant columns
  fold_cols <-  grep(":log2FoldChange", colnames(linda.output_reduced), value= TRUE)
  
  pval_cols <- grep(":padj", colnames(linda.output_reduced), value= TRUE)
  
  # Create a data frame with estimates and p-values
  plot_df <- linda.output_reduced %>%
    dplyr::select("SeqID",all_of(fold_cols), dplyr::all_of(pval_cols)) %>%
    dplyr::mutate(across(ends_with(":padj"), ~ ifelse(. < 0.001, "***", ifelse(. < 0.01, "**", ifelse(. < 0.05, "*", "")))))
  
  plot_df_melt_fold <- melt(plot_df, id.vars = "SeqID", measure.vars = fold_cols, variable.name = "Variable", value.name = "log2FoldChange")
  plot_df_melt_pval <- melt(plot_df, id.vars = "SeqID", measure.vars = pval_cols, variable.name = "Variable", value.name = "Significance")
  
  # Remove ":Est" and ":pval" from variable names for easier merging
  plot_df_melt_fold$Variable <- gsub(":log2FoldChange", "", plot_df_melt_fold$Variable)
  plot_df_melt_pval$Variable <- gsub(":padj", "", plot_df_melt_pval$Variable)
  
  # Combine estimates and significance stars into one dataframe
  plot_df_combined <- base::merge(plot_df_melt_fold, plot_df_melt_pval, by = c("SeqID","Variable"))
  
  plot_df_combined$fold <- as.numeric(as.character(plot_df_combined$log2FoldChange))
  
  plot_df_combined <- plot_df_combined %>%
    mutate(log2FoldChange = ifelse(is.finite(log2FoldChange), log2FoldChange, NA))
  
  # Check for NAs and remove them if necessary
  #plot_df_combined <- na.omit(plot_df_combined)
  
  # Pivot the combined data for heatmap
  taxa_tab %<>% column_to_rownames("SeqID")
  heatmap_data <- maditr::dcast(plot_df_combined, SeqID ~ Variable, value.var = "log2FoldChange")
  if(is_dna_sequence(heatmap_data$SeqID[1])){
    heatmap_taxa <- taxa_tab[heatmap_data$SeqID,]
    
    # genus assigning
    genus <- heatmap_taxa$Genus
    
    # family assigning
    where_unassigned <- grep("unassigned",genus)
    where_uncultured <- grep("uncultured",genus)
    where_unassigned <- c(where_unassigned,where_uncultured)
    genus[where_unassigned] <- paste0("f__",heatmap_taxa$Family[where_unassigned],";g__",genus[where_unassigned])
    
    # order assigning
    where_f_unassigned <-  grep("f__unassigned",genus)
    where_f_uncultured <- grep("f__uncultured",genus)
    where_f_unassigned <- c(where_f_unassigned,where_f_uncultured)
    genus[where_f_unassigned] <- paste0("o__",heatmap_taxa$Order[where_f_unassigned],";",genus[where_f_unassigned])
    
    duplicated <- genus[duplicated(genus)]
    
    for (duplic in unique(duplicated)){
      where_duplic <- which(genus==duplic)
      genus[where_duplic] <- paste(duplic,1:length(where_duplic))
    }
  } else genus <- heatmap_data$SeqID
  
  
  heatmap_data$SeqID <- genus
  #heatmap_data$SeqID[duplicated(heatmap_data$SeqID)] <- paste(duplicated,1:length(duplicated))
  heatmap_data %<>% column_to_rownames(var="SeqID")
  
  # Prepare annotation for significance stars
  annotation <- dcast(plot_df_combined, SeqID ~ Variable, value.var = "Significance")
  annotation %<>% mutate(SeqID=rownames(heatmap_data)) %>% column_to_rownames("SeqID") 
  annotation[is.na(annotation)] <- ""
  # Plot the heatmap
  p <- as.ggplot(pheatmap(heatmap_data , 
                          cluster_rows = FALSE, 
                          cluster_cols = FALSE, 
                          display_numbers = annotation, 
                          number_color = "black",
                          fontsize_number = 8,
                          main = "LinDA's log2FoldChange and Significance",silent=TRUE,
                          angle_col=0))
  
  return(p)
}

dot_heatmap_linda <- function(raw_linda, uni_list,
                              taxa_table, group=NULL){
  if (class(raw_linda)=="data.frame") {
    raw_linda <- raw_linda[, c("SeqID","Taxonomy","padj", "log2FoldChange","MEDIAN_clr_ALL")]
    raw_linda <- raw_linda[,c(TRUE,!grepl("ASV",colnames(raw_linda)[-1]))]
    raw_linda <- raw_linda[,c(TRUE,!grepl("Taxonomy",colnames(raw_linda)[-2]))]
    colnames(raw_linda) <- c("SeqID","Taxonomy",paste0(group,c(":padj",":log2FoldChange",":MEDIAN_clr_ALL")))
  }
  else {
    if (is.null(group)){
      if (TRUE %in% grepl("Group",names(raw_linda))) wanted_list <- raw_linda[grepl("vs Group",names(raw_linda)) & !grepl(":",names(raw_linda))]
      else wanted_list <- raw_linda
    } else wanted_list <- raw_linda[group]
    
    # Assuming wanted_list and uni_list are lists of dataframes
    raw_linda <- Map(function(wanted_df, uni_df) {
      merge(wanted_df, uni_df[, c("SeqID", "MEDIAN_clr_ALL")], by = "SeqID", all.x = TRUE)
    }, wanted_list, uni_list)
    
    
    # raw_linda <- lapply(wanted_list, function(df) {
    #   df <- merge(df,uni_df[,c("SeqID","MEDIAN_clr_ALL")],by="SeqID",all.x = TRUE)
    #   df[, c("SeqID","Taxonomy","padj", "log2FoldChange","MEDIAN_clr_ALL")]
    #   return(df)
    #   })
    #raw_linda <- lapply(wanted_list, function(df) df[, c("ASV","Taxonomy","padj", "log2FoldChange","MEDIAN ALL","Taxonomy")])
    
    names(raw_linda) <- gsub("(terminal_ileum)|(colon)|(ASV)|(genus)","",names(raw_linda))
    names(raw_linda) <- gsub("^ +","",names(raw_linda))
    
    groups <- names(raw_linda)
    groups <- rep(groups, each = 3)
    
    raw_linda <- 
      suppressWarnings(Reduce(function(x, y) merge(x, y, by = "SeqID", all = TRUE), raw_linda))
    
    if (is_dna_sequence(raw_linda$SeqID[1])){
      for (i in 1:nrow(raw_linda)){
        where_taxonomy <- grep("Taxonomy", colnames(raw_linda))
        where_nonnan <- grep(TRUE,(!is.na(raw_linda[i,where_taxonomy])))[1]
        nonnan_taxonomy <- raw_linda[i,where_taxonomy[where_nonnan]]
        nonnan_taxonomy_genus <- regmatches(nonnan_taxonomy, regexpr("g__.+;s__", nonnan_taxonomy))
        nonnan_taxonomy_genus <- substring(nonnan_taxonomy_genus,4,nchar(nonnan_taxonomy_genus)-4)
        if ((nonnan_taxonomy_genus == "unassigned")|(nonnan_taxonomy_genus == "uncultured")){
          nonnan_taxonomy_genus <- regmatches(nonnan_taxonomy, regexpr("f__.+;s__", nonnan_taxonomy))
          nonnan_taxonomy_genus <- substring(nonnan_taxonomy_genus,1,nchar(nonnan_taxonomy_genus)-4)
        }
        raw_linda[i,"SeqID"] <- nonnan_taxonomy_genus
      }
    }
    
    asvs <- raw_linda$SeqID
    duplicated <- asvs[duplicated(asvs)]
    
    for (duplic in unique(duplicated)){
      where_duplic <- which(asvs==duplic)
      asvs[where_duplic] <- paste(duplic,1:length(where_duplic))
    }
    raw_linda$SeqID <- asvs 
    
    raw_linda <- raw_linda[,c(TRUE,!grepl("SeqID",colnames(raw_linda)[-1]))]
    raw_linda <- raw_linda[,c(TRUE,!grepl("Taxonomy",colnames(raw_linda)[-2]))]
    raw_linda <- raw_linda[,!grepl("p_value",colnames(raw_linda))]
    
    colnames(raw_linda) <- c("SeqID","Taxonomy",paste0(groups,c(":log2FoldChange",":padj",":MEDIAN_clr_ALL")))
  }
  #raw_linda %<>% dplyr::bind_cols(.name_repair = "unique")
  
  #taxa_table <- taxa_table[taxa_table$SeqID %in% raw_linda$SeqID,]
  #raw_linda$Taxonomy=create_asv_taxa_table(raw_linda,taxa_table)$SeqID
  
  ##########################################################################################################
  #raw_linda[,grepl(":padj",colnames(raw_linda))][is.na(raw_linda[,grepl(":padj",colnames(raw_linda))])] <- Inf
  #raw_linda[is.na(raw_linda)] <- 0
  
  or_logic <- rep(FALSE,dim(raw_linda)[1])
  for (col in grep(":padj",colnames(raw_linda))){
    or_logic <- (or_logic) | ifelse(is.na(raw_linda[,col] < 0.1),
                                    FALSE,
                                    (raw_linda[,col] < 0.1))
  }
  
  raw_linda_reduced <- raw_linda[or_logic,] # %>%
  #  rownames_to_column("SeqID")
  
  # Extract relevant columns
  fold_cols <-  grep(":log2FoldChange", colnames(raw_linda_reduced), value= TRUE)
  pval_cols <- grep(":padj", colnames(raw_linda_reduced), value= TRUE)
  median_cols <- grep(":MEDIAN_clr_ALL", colnames(raw_linda_reduced), value= TRUE)
  effect_columns <- grep("effect:padj", colnames(raw_linda_reduced))
  
  # Create a data frame with estimates and p-values
  plot_df <- raw_linda_reduced %>%
    dplyr::select("SeqID","Taxonomy",all_of(fold_cols), dplyr::all_of(pval_cols),dplyr::all_of(median_cols)) %>%
    dplyr::mutate(across(ends_with("effect"), ~ ifelse(. == 999, "\u26AB", ""))) %>%
    dplyr::mutate(across(ends_with(":padj"), ~ ifelse(.==999, "+",ifelse(. < 0.001, "***", ifelse(. < 0.01, "**", ifelse(. < 0.05, "*", ""))))) )


  if (is_dna_sequence(plot_df$SeqID[1])){
    # genus assigning
    #regmatches(plot_df$Taxonomy, regexpr("g__.+;", plot_df$Taxonomy)) 
    genus <- regmatches(plot_df$Taxonomy, regexpr("g__.+;", plot_df$Taxonomy)) 
    genus <- substring(genus,4,nchar(genus)-1)
    
    # family assigning
    where_unassigned <- grep("unassigned",genus)
    where_uncultured <- grep("uncultured",genus)
    where_unassigned <- c(where_unassigned,where_uncultured)
    family <- regmatches(plot_df$Taxonomy[where_unassigned], regexpr("f__.+;", plot_df$Taxonomy[where_unassigned])) 
    family <- substring(family,1,nchar(family)-1)
    genus[where_unassigned] <- family
    
    # order assigning
    where_f_unassigned <-  grep("f__unassigned",genus)
    where_f_uncultured <- grep("f__uncultured",genus)
    where_f_unassigned <- c(where_f_unassigned,where_f_uncultured)
    order <- regmatches(plot_df$Taxonomy[where_f_unassigned], regexpr("o__.+;", plot_df$Taxonomy[where_f_unassigned])) 
    order <- substring(order,1,nchar(order)-1)
    genus[where_f_unassigned] <- order
    
    duplicated <- genus[duplicated(genus)]
    
    for (duplic in unique(duplicated)){
      where_duplic <- which(genus==duplic)
      genus[where_duplic] <- paste(duplic,1:length(where_duplic))
    }
    
    plot_df$SeqID <- genus
  }
  
  # ZORADENIE PODLA LOGFOLDCHANGE
  values <- plot_df[,c(1,3)]
  values[,2][is.na(values[,2])] <- plot_df[is.na(values[,2]),4]
  values[,2][is.na(values[,2])] <- plot_df[is.na(values[,2]),5]
  ordering <- order(values[,2])
  
  plot_df <- plot_df[ordering,]
  rownames(plot_df) <- NULL
  
  plot_df_melt_fold <- melt(plot_df, id.vars = "SeqID", measure.vars = fold_cols, variable.name = "Variable", value.name = "log2FoldChange")
  plot_df_melt_pval <- melt(plot_df, id.vars = "SeqID", measure.vars = pval_cols, variable.name = "Variable", value.name = "Significance")
  plot_df_melt_median <- melt(plot_df, id.vars = "SeqID", measure.vars = median_cols, variable.name = "Variable", value.name = "Median")
  
  # Remove ":Est" and ":pval" from variable names for easier merging
  plot_df_melt_fold$Variable <- gsub(":log2FoldChange", "", plot_df_melt_fold$Variable)
  plot_df_melt_pval$Variable <- gsub(":padj", "", plot_df_melt_pval$Variable)
  plot_df_melt_median$Variable <- gsub(":MEDIAN_clr_ALL", "", plot_df_melt_median$Variable)
  
  # Combine estimates and significance stars into one dataframe
  plot_df_combined <- base::merge(plot_df_melt_fold, plot_df_melt_pval, by = c("SeqID","Variable"))
  plot_df_combined <- base::merge(plot_df_combined, plot_df_melt_median, by = c("SeqID","Variable"))
  
  # Check for NAs and remove them if necessary
  #plot_df_combined <- na.omit(plot_df_combined)
  
  plot_df_combined %<>% mutate(median_clr = Median)
  plot_df_combined$SeqID <- factor(plot_df_combined$SeqID,levels = plot_df$SeqID)
  min_clr <- min(plot_df_combined$median_clr, na.rm = TRUE)
  cat("min_clr",min_clr,"\n")
  max_clr <- max(plot_df_combined$median_clr, na.rm = TRUE)
  cat("max_clr",max_clr,"\n")
  
  min_log <- min(plot_df_combined$log2FoldChange, na.rm = TRUE)
  cat("min_log",min_log,"\n")
  max_log <- max(plot_df_combined$log2FoldChange, na.rm = TRUE)
  cat("max_log",max_log,"\n")
  
  plot_df_combined$Variable <- gsub("pre_ltx","pre_LTx",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("post_ltx","post_LTx",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("healthy","HC",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("HC vs pre_LTx","pre_LTx vs HC",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("HC vs post_LTx","post_LTx vs HC",plot_df_combined$Variable)
  
  plot_df_combined$Variable <- gsub("HC vs rPSC","rPSC vs HC",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("HC vs non-rPSC","non-rPSC vs HC",plot_df_combined$Variable)
  plot_df_combined$Variable <- gsub("non-rPSC vs rPSC","rPSC vs non-rPSC",plot_df_combined$Variable)
  
  plot_df_combined$Variable <- gsub("pre_LTx vs post_LTx","post_LTx vs pre_LTx",plot_df_combined$Variable)

  plot_df_combined$Variable <- gsub("  ","",plot_df_combined$Variable)
  
  if (TRUE %in% (grepl("post_LTx",plot_df_combined$Variable))){
    plot_df_combined$Variable <- factor(plot_df_combined$Variable,
                                        levels=c("pre_LTx vs HC","post_LTx vs HC","post_LTx vs pre_LTx"))
  } else {
    plot_df_combined$Variable <- factor(plot_df_combined$Variable,
                                        levels=c("non-rPSC vs HC","rPSC vs HC","rPSC vs non-rPSC"))
  }
  
  p <- ggplot(plot_df_combined, aes(x = Variable, y = SeqID)) +
    geom_point(aes(size = median_clr, fill = log2FoldChange, color=log2FoldChange), shape = 21) +
    geom_text(aes(label = Significance,vjust=0.6), color = "black", size = 2) +  # Add asterisks for significant points
    #geom_text(aes(label = effect,vjust=0.6), color = "black", size = 2,) +
    #scale_size_continuous(name = "Median clr") +
    scale_size_continuous(name = "Median clr",range = c(1, 5), 
                          limits = c(-5, 7.1),
                          breaks = c(-5, -2.5, 0, 2.5, 5)) +
    scale_fill_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-5.5, 5.5)) +
    scale_color_gradient2(name = "log2FoldChange", low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(-5.5, 5.5)) +
    #scale_color_manual(values = "black")
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
    #labs(title = "LinDA's log2FoldChange, Significance and median clr value") 
  
  return(p)
  
}

volcano_plot_linda <- function(linda.output,group,taxa_table,cutoff.pval=0.05, cutoff.lfc=1){
  
  output <- linda.output[[group]]
  #output <- output[taxa_table$SeqID,]
  
  lfc <- output$log2FoldChange
  padj <- output$padj
  
  called <- output[,"padj"] <= cutoff.pval
  all.p <- output[,"padj"]
  all.col <- rgb(0, 0, 0, 0.2)
  
  data_df <- data.frame(x=lfc,
                        y=-1*log10(all.p),
                        name=rownames(output))
  
  data_called <- data_df[called,]
  taxa_table <- taxa_table %>% column_to_rownames("SeqID")
  if ("Genus" %in% colnames(taxa_table)){
    data_called$name <- taxa_table[data_called$name,"Genus"]
  } else {
    message("Using Phylum for naming")
    data_called$name <- taxa_table[data_called$name,"Phylum"]
  }
  
  
  maximum <- max(c(abs(min(data_df$x)), abs(max(data_df$x))))
  p <- ggplot(data=data_df, aes(x=x,y=y)) +
    geom_point(colour=all.col, shape=19, size=3) + theme_bw() + 
    xlab("Log"[2]~"Fold Change ") + ylab("-1 * Median Log"[10]~" q value") +
    geom_point(data=data_called, aes(x=x,y=y), colour="red", size=3)+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15),
          axis.title=element_text(size=10)) + 
    geom_vline(xintercept=1.5, linetype=2, colour="grey")+ 
    geom_vline(xintercept=-1.5, linetype=2, colour="grey")+ 
    coord_cartesian(xlim = c(-maximum,maximum), clip="on") +
    geom_hline(yintercept=-1*log10(0.05),linetype=2, colour="grey") +
    geom_text_repel(data=data_called, aes(x = x, y = y, label = name),size=2,max.overlaps = 20)
  
  return(p)
}

volcano_plot_ancom <- function(ancom_output,taxa_table,cutoff.pval=0.05, cutoff.lfc=1){
  
  output <- ancom_output
  
  lfc <- output$LFC
  padj <- output$p.adjusted
  
  called <- output[,"p.adjusted"] <= cutoff.pval
  all.p <- output[,"p.adjusted"]
  all.col <- rgb(0, 0, 0, 0.2)
  
  data_df <- data.frame(x=lfc,
                        y=-1*log10(all.p),
                        name=output$SeqID)
  
  data_called <- data_df[called,]
  taxa_table <- taxa_table %>% column_to_rownames("SeqID")
  data_called$name <- taxa_table[data_called$name,"Genus"]
  
  maximum <- max(c(abs(min(data_df$x)), abs(max(data_df$x))))
  p <- ggplot(data=data_df, aes(x=x,y=y)) +
    geom_point(colour=all.col, shape=19, size=3) + theme_bw() + 
    xlab("Log"[2]~"Fold Change ") + ylab("-1 * Median Log"[10]~" q value") +
    geom_point(data=data_called, aes(x=x,y=y), colour="red", size=3)+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15),
          axis.title=element_text(size=10)) + 
    geom_vline(xintercept=1.5, linetype=2, colour="grey")+ 
    geom_vline(xintercept=-1.5, linetype=2, colour="grey")+ 
    coord_cartesian(xlim = c(-maximum,maximum), clip="on") +
    geom_hline(yintercept=-1*log10(0.05),linetype=2, colour="grey") +
    geom_text_repel(data=data_called, aes(x = x, y = y, 
                                          label = name),size=2,max.overlaps = 20)
  
  return(p)
}


volcano_plot_maaslin <- function(maaslin_output,taxa_table,cutoff.pval=0.05, cutoff.lfc=1,variable="Group"){
  
  # group effect
  output <- maaslin_output$results
  output <- output[output$metadata==variable,c(1,4,5,6,8)]
  
  lfc <- output$coef
  padj <- output$qval
  
  called <- output[,"qval"] <= cutoff.pval
  all.p <- output[,"qval"]
  all.col <- rgb(0, 0, 0, 0.2)
  
  data_df <- data.frame(x=lfc,
                        y=-1*log10(all.p),
                        name=output$feature)
  
  data_df[is.na(data_df$x),] <- 0
  
  data_called <- data_df[called,]
  taxa_table <- taxa_table %>% column_to_rownames("SeqID")
  taxonomic_level <- colnames(taxa_table)[ncol(taxa_table)]
  if (taxonomic_level=="Species") taxonomic_level = "Genus"
  data_called$name <- taxa_table[data_called$name,taxonomic_level]
  
  maximum <- max(c(abs(min(data_df$x)), abs(max(data_df$x))))
  p <- ggplot(data=data_df, aes(x=x,y=y)) +
    geom_point(colour=all.col, shape=19, size=3) + theme_bw() + 
    xlab("Maaslin2 coefficient") + ylab("-1 * Median Log"[10]~" q value") +
    geom_point(data=data_called, aes(x=x,y=y), colour="red", size=3)+
    theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15),
          axis.title=element_text(size=10)) + 
    geom_vline(xintercept=1.5, linetype=2, colour="grey")+ 
    geom_vline(xintercept=-1.5, linetype=2, colour="grey")+ 
    coord_cartesian(xlim = c(-maximum,maximum), clip="on") +
    geom_hline(yintercept=-1*log10(0.05),linetype=2, colour="grey") +
    geom_text_repel(data=data_called, aes(x = x, y = y, 
                                          label = name),size=2,max.overlaps = 20)

  return(p)

}



roc_curve <- function(enet_model, group){
  ggroc_data <- ggroc(enet_model$kfold_rocobjs)$data
  roc_c <- ggplot(data=ggroc_data) + 
    geom_line(aes(x=`1-specificity`, y=sensitivity, by=name, color="red",alpha=0.9)) +
    #geom_smooth(aes(x=`1-specificity`, y=sensitivity),
    #            linewidth=1,
    #            color = 'red',
    #            se=FALSE,
    #            method = 'loess') +
    theme_minimal() + 
    theme(legend.position = "none") + 
    ggtitle(paste0('ROC Curve: ',group[1],' vs ',group[2],' (AUC = ', round(enet_model$model_summary$auc_optimism_corrected,2), ')'))  
    return(roc_c)
}

roc_curve_all <- function(objects){
  p <- ggplot()
  colors <- c(
    "red", "blue","#2B9D2B", "#D90368",
    "#984ea3", "#FAF33E","#DF75BE", "grey",
    "#17BACB", "#66c2a5","#A5BE00", "#000000",
    "#a65629")
  
  colors <- c("#4169E1","#984ea3","#008080",
              "#FF6347","#FFD700","#D90368")
  
  
  #colors <- c(
  #  "#f6bcb7","#d8d082","#88db9c",
  #  "#8cdee0","#b5ccfe","#f6b2f0")
  

  names(colors) <- names(objects)
  
  for (i in 1:length(objects)){
    auc <- objects[[i]]$auc
    my_df <- data.frame(sensitivity=objects[[i]]$sensitivities,
                        `1-specificity`=1-objects[[i]]$specificities,
                        check.names = FALSE)
  
      
    my_color <- names(colors)[i]
    p <- p + 
    geom_line(data=my_df,aes(x=`1-specificity`, y=sensitivity, group=name,
                                  color=!!my_color), linewidth=1.5,alpha=1) 


  }
  p <- p + theme_minimal() + 
    scale_color_manual(values = colors) + #+ guides(colour = "none") + 
  theme(legend.title = element_blank())
  return(p)
}


roc_curve_all_custom <- function(objects,Q,model_name,legend=TRUE){
  print(names(objects))
  p <- ggplot()
  #colors <- c(
  #  "red", "blue","#2B9D2B", "#D90368",
  #  "#984ea3", "#FAF33E","#DF75BE", "grey",
  #  "#17BACB", "#66c2a5","#A5BE00", "#000000",
  #  "#a65629")
  
  #colors <- c(
  #  "#708090", "#17BACB","#FFA500", "#7B68EE",
  #  "#FF6347", "#F0E68C","#708090")
  
  #colors <- c("#FF6347","#556B2F","#7B68EE","#17BACB","#708090","#FFA500")
  
  if (grepl("^Q1", Q)) colors <- c("#FF6347","#708090","#17BACB")
  if (grepl("^Q2",Q)) {
    #colors <- c("#FF7F50","#FFD700","#4169E1","#008080") 
    colors <- c("#008080","#FFD700","#4169E1","#FF7F50")
  }
  if (grepl("^Q3",Q)) {
    #\colors <- c("#FF7F50","#FFD700","#4169E1","#008080") 
    colors <- c("#000000")
  }
  
  names(colors) <- names(objects)
  
  for (i in 1:length(objects)){
    loaded_name <- load(file.path("../intermediate_files","models",Q,names(objects)[i],paste0(model_name,".RData")))
    my_model <- get(loaded_name)
    auc_optimism_corrected <- my_model$model_summary$auc_optimism_corrected
                  
    aucs <- sapply(objects[[i]],function(df) df$auc)
    q1_auc <- quantile(aucs, probs = 0.025)
    q3_auc <- quantile(aucs, probs = 0.975)
    
    differences <- abs(aucs - q1_auc)
    min_auc <- which.min(differences)
    
    differences <- abs(aucs - q3_auc)
    max_auc <- which.min(differences)
    
    differences <- abs(aucs - auc_optimism_corrected)
    true_auc <- which.min(differences)
    
    #a <- objects[[i]][[max_auc]]
    #b <- objects[[i]][[min_auc]]
    c <- objects[[i]][[true_auc]]
    ggroc_data_c <- ggroc(c)$data
    #y3_interp <- approx(ggroc_data_c$`1-specificity`, ggroc_data_c$sensitivity, xout = x_common)$y
    
    ############################################
    # Create a common set of x-axis points
    x_common <- seq(0, 1, length.out = 100)
    
    my_df <- data.frame(x=x_common)
    for (j in 1:length(aucs)){
      ggroc_data <- ggroc(objects[[i]][[j]])$data
      ggroc_data$`1-specificity` <- ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/100000000000000000000
      ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
      y_aprox <- approx(ggroc_data$`1-specificity`, ggroc_data$sensitivity, xout = x_common)$y
      my_df <- cbind(my_df,y_aprox)
    }
    my_df <- my_df %>% column_to_rownames("x")
    q1_y <- c()
    q2_y <- c()
    q3_y <- c()
    for (j in 1:length(x_common)){
      q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), probs = 0.025,na.rm=TRUE))
      q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
      q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), probs = 0.975,na.rm=TRUE))
    }
    df <- data.frame(x = x_common, 
                     y1 = q1_y, 
                     y2 = q3_y, 
                     y3 = q2_y)
    df[is.na(df)] <- 0
    ##########################################################
    # Interpolate y-values for both lines
    #y1_interp <- approx(ggroc_data_a$`1-specificity`, ggroc_data_a$sensitivity, xout = x_common)$y
    #y2_interp <- approx(ggroc_data_b$`1-specificity`, ggroc_data_b$sensitivity, xout = x_common)$y
    #y3_interp <- approx(ggroc_data_c$`1-specificity`, ggroc_data_c$sensitivity, xout = x_common)$y
    
    #df <- data.frame(x = x_common, 
    #                 y1 = y1_interp, 
    #                 y2 = y2_interp, 
    #                 y3 = y3_interp)
    
    my_color <- colors[i]
    
    #p <- p + 
    #geom_line(data=ggroc_data,aes(x=`1-specificity`, y=sensitivity, group=name,
    #                              color=!!my_color), linewidth=2,alpha=0.7) +
    #geom_ribbon(data = bounds, aes(x, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.4) 
    #geom_smooth(data=ggroc_data,aes(x=`1-specificity`, y=sensitivity,color = !!my_color),
    #            linewidth=1,
    #            se=FALSE,
    #            method = 'loess')
    p <- p + 
      #  geom_line(data=df,aes(x=x,y = y1,color=!!my_color)) +
      geom_ribbon(data=df,aes(x =x,ymin = pmin(y1, y2), ymax = pmax(y1, y2)), fill=my_color, alpha = 0.5,show.legend = TRUE) +
      #geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1.5) +
      theme_classic() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            axis.title.x = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.y = element_text(size=5)) + 
      ylab("Sensitivity") + xlab("1-specificity")
  }
  for (i in 1:length(objects)){
    loaded_name <- load(file.path("../intermediate_files","models",Q,names(objects)[i],paste0(model_name,".RData")))
    my_model <- get(loaded_name)
    auc_optimism_corrected <- my_model$model_summary$auc_optimism_corrected
    
    aucs <- sapply(objects[[i]],function(df) df$auc)
    q1_auc <- quantile(aucs, probs = 0.025)
    q3_auc <- quantile(aucs, probs = 0.975)
    
    differences <- abs(aucs - q1_auc)
    min_auc <- which.min(differences)
    
    differences <- abs(aucs - q3_auc)
    max_auc <- which.min(differences)
    
    differences <- abs(aucs - auc_optimism_corrected)
    true_auc <- which.min(differences)
    
    #a <- objects[[i]][[max_auc]]
    #b <- objects[[i]][[min_auc]]
    c <- objects[[i]][[true_auc]]
    ggroc_data_c <- ggroc(c)$data
    y3_interp <- approx(ggroc_data_c$`1-specificity`, ggroc_data_c$sensitivity, xout = x_common)$y
    
    ############################################
    # Create a common set of x-axis points
    x_common <- seq(0, 1, length.out = 100)
    
    my_df <- data.frame(x=x_common)
    for (j in 1:length(aucs)){
      ggroc_data <- ggroc(objects[[i]][[j]])$data
      ggroc_data$`1-specificity` <- ggroc_data$`1-specificity` + (1:nrow(ggroc_data))/100000000000000000000
      ggroc_data$`1-specificity`[nrow(ggroc_data)] <- 1
      y_aprox <- approx(ggroc_data$`1-specificity`, ggroc_data$sensitivity, xout = x_common)$y
      my_df <- cbind(my_df,y_aprox)
    }
    my_df <- my_df %>% column_to_rownames("x")
    q1_y <- c()
    q2_y <- c()
    q3_y <- c()
    for (j in 1:length(x_common)){
      q1_y <- c(q1_y,quantile(as.numeric(my_df[j,]), probs = 0.025,na.rm=TRUE))
      q2_y <- c(q2_y,mean(as.numeric(my_df[j,])))
      q3_y <- c(q3_y,quantile(as.numeric(my_df[j,]), probs = 0.975,na.rm=TRUE))
    }
    df <- data.frame(x = x_common, 
                     y1 = q1_y, 
                     y2 = q3_y, 
                     y3 = q2_y)
    df[is.na(df)] <- 0
    ##########################################################
    # Interpolate y-values for both lines
    #y1_interp <- approx(ggroc_data_a$`1-specificity`, ggroc_data_a$sensitivity, xout = x_common)$y
    #y2_interp <- approx(ggroc_data_b$`1-specificity`, ggroc_data_b$sensitivity, xout = x_common)$y
    #y3_interp <- approx(ggroc_data_c$`1-specificity`, ggroc_data_c$sensitivity, xout = x_common)$y
    
    #df <- data.frame(x = x_common, 
    #                 y1 = y1_interp, 
    #                 y2 = y2_interp, 
    #                 y3 = y3_interp)
    
    my_color <- colors[i]
    
    #p <- p + 
    #geom_line(data=ggroc_data,aes(x=`1-specificity`, y=sensitivity, group=name,
    #                              color=!!my_color), linewidth=2,alpha=0.7) +
    #geom_ribbon(data = bounds, aes(x, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.4) 
    #geom_smooth(data=ggroc_data,aes(x=`1-specificity`, y=sensitivity,color = !!my_color),
    #            linewidth=1,
    #            se=FALSE,
    #            method = 'loess')
    p <- p + 
      #  geom_line(data=df,aes(x=x,y = y1,color=!!my_color)) +
      #geom_ribbon(data=df,aes(x =x,ymin = pmin(y1, y2), ymax = pmax(y1, y2)), fill=my_color, alpha = 0.5,show.legend = TRUE) +
      geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1) +
      theme_classic() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            axis.title.x = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.y = element_text(size=5)) + 
      ylab("Sensitivity") + xlab("1-specificity")
  }
  
  if (legend){
    if (TRUE %in% grepl("post_ltx",names(objects))){
      legend_data <- data.frame(
        color = gsub("(terminal_ileum)|(colon)|(Genus)|(ASV)","",names(objects)),
        value = c(3,1,2),
        color_hex <- colors
      )
    } else if (TRUE %in% grepl("rPSC",names(objects))){
      legend_data <- data.frame(
        color = gsub("(terminal_ileum)|(colon)|(Genus)|(ASV)","",names(objects)),
        value = c(2,3,4,1),
        color_hex <- colors)
    } else{
      legend_data <- data.frame(
        color = gsub("(terminal_ileum)|(colon)|(Genus)|(ASV)","",names(objects)),
        value = c(1),
        color_hex <- colors)
    }
    
    
    legend_data$color <- gsub("pre_ltx","pre_LTx",legend_data$color)
    legend_data$color <- gsub("post_ltx","post_LTx",legend_data$color)
    legend_data$color <- gsub("healthy","HCs",legend_data$color)
    
    # Custom legend plot
    legend_plot <- ggplot(legend_data, aes(x = 0.5, y = value, fill = color_hex)) +
      geom_tile(width = 0.3, height = 0.8) +
      geom_text(size=4,aes(label = color), hjust = 0, nudge_x = 0.2, color = "black") +
      scale_fill_identity() +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) + xlim(0, 5)  # Adjust spacing between tiles and text
    
    p <- ggarrange(p, legend_plot, ncol=1,heights = c(1,0.2))  
  }
  
  return(p)
}


## ANCOM functions ----

ancom_analysis <- function(asv_table,taxa_table,metadata, fix_formula="",random_formula=NULL,
                           prv_cut=0.0,lib_cut=0,s0_perc=0){
  if (FALSE %in% ("SeqID" %in% colnames(asv_table))) asv_table %<>% rownames_to_column("SeqID")
  if (FALSE %in% ("SampleID" %in% colnames(metadata))) metadata %<>% rownames_to_column("SampleID")
  phyloseq_for_ancom <- construct_phyloseq(asv_table,taxa_table,metadata)
  output <- ANCOMBC::ancombc2(data = phyloseq_for_ancom,
                              fix_formula = fix_formula, rand_formula = NULL,
                              p_adj_method = "BH", pseudo_sens = TRUE,
                              prv_cut = prv_cut, lib_cut = lib_cut, s0_perc = s0_perc,
                              group = "Group",
                              alpha = 0.05, verbose = TRUE,
                              mdfdr_control = list(fwer_ctrl_method = "BH", B = 100))
  #tryCatch(
  #  output <- ANCOMBC::ancombc2(data = phyloseq_for_ancom,
  #                    fix_formula = fix_formula, rand_formula = NULL,
  #                    p_adj_method = "BH", pseudo_sens = TRUE,
  #                    prv_cut = prv_cut, lib_cut = lib_cut, s0_perc = s0_perc,
  #                    group = "Group",
  #                    alpha = 0.05, verbose = TRUE,
  #                    mdfdr_control = list(fwer_ctrl_method = "BH", B = 100)),
  #  error = function(cnd) {
  #    bad_taxa <- strsplit(conditionMessage(cnd), "\n")[[c(1L, 2L)]]
  #    bad_taxa <- .subset2(strsplit(bad_taxa, ", "), 1L)
  #    return(bad_taxa)
  #  }
  #)
  #return(bad_taxa)

  
  result <- output$res[,c(TRUE,grepl("Group",colnames(output$res))[-1])]
  colnames(result) <- c("SeqID","LFC","SE","w","p.value","p.adjusted",
                        "significance","passed_ss")
  return(result)
}

raw_ancom.df <- function(ancom_output,uni_data,uni_taxa){
  
  if (is_dna_sequence(rownames(uni_data)[1])) {
  raw_ancom_result <- data.frame(
    ASV=ancom_output$SeqID,
    Taxonomy=create_asv_taxa_table(uni_data[ancom_output$SeqID,] %>% rownames_to_column("SeqID"),uni_taxa)$SeqID,
    log2FoldChange=ancom_output$LFC,
    p_value=ancom_output$p.value,
    padj=ancom_output$p.adjusted)
  } else {
    raw_ancom_result <- data.frame(
      ASV=ancom_output$SeqID,
      log2FoldChange=ancom_output$LFC,
      p_value=ancom_output$p.value,
      padj=ancom_output$p.adjusted)
  }
  
  return(raw_ancom_result)
}

raw_maaslin.df <- function(fit_data,uni_data,uni_taxa){
  fit_data <- fit_data$results
  fit_data <- fit_data[fit_data$metadata=="Group",c(1,4,5,6,8)]

  if (is_dna_sequence(rownames(uni_data)[1])) {
    raw_maaslin_result <- data.frame(
      ASV=fit_data$feature,
      Taxonomy=create_asv_taxa_table(uni_data[fit_data$feature,] %>% rownames_to_column("SeqID"),uni_taxa)$SeqID,
      coef=fit_data$coef,
      p_value=fit_data$pval,
      padj=fit_data$qval)
  } else {
    raw_maaslin_result <- data.frame(
      ASV=fit_data$feature,
      coef=fit_data$coef,
      p_value=fit_data$pval,
      padj=fit_data$qval)
  }
  
  return(raw_maaslin_result)
}


## linDA functions ----

rawlinda.df <- function(linda.output,group,uni_data,uni_taxa){
  if (is_dna_sequence(rownames(uni_data)[1])) {
  raw_linda_result <- data.frame(
    SeqID=rownames(linda.output[[group]]),
    Taxonomy=create_asv_taxa_table(uni_data[rownames(linda.output[[group]]),] %>% 
                                     rownames_to_column("SeqID"),uni_taxa)$SeqID,
    log2FoldChange=linda.output[[group]]$log2FoldChange,
    p_value=linda.output[[group]]$pvalue,
    padj=linda.output[[group]]$padj)
  } else {
    raw_linda_result <- data.frame(
      SeqID=rownames(linda.output[[group]]),
      Taxonomy=create_asv_taxa_table(uni_data[rownames(linda.output[[group]]),] %>% 
                                       rownames_to_column("SeqID"),uni_taxa)$SeqID,
      log2FoldChange=linda.output[[group]]$log2FoldChange,
      p_value=linda.output[[group]]$pvalue,
      padj=linda.output[[group]]$padj)
  }
  return(raw_linda_result)
}

linda.df <- function(linda.output,group,filt_ileum_uni_data,filt_ileum_uni_taxa){
  filt_ileum_uni_taxa <- filt_ileum_uni_taxa %>% column_to_rownames("SeqID")
  diff_asvs <- rownames(linda.output[[group]])[linda.output[[group]]$padj < 0.05]
  
  # save the results
  linda_result <- data.frame(
    ASVs=diff_asvs,
    Taxonomy=create_asv_taxa_table(
      filt_ileum_uni_data[diff_asvs,] %>% rownames_to_column("SeqID"),
      filt_ileum_uni_taxa[diff_asvs,] %>% rownames_to_column("SeqID"))$SeqID,
    log2FoldChange=linda.output[[group]][diff_asvs,"log2FoldChange"],
    pvalue=linda.output[[group]][diff_asvs,"pvalue"],
    padj=linda.output[[group]][diff_asvs,"padj"])
  
  return(linda_result)
}



linda_renaming <- function(linda_data, group){
  # This function renames the results of linDA (linda() function $output) based on the chosen intercept 
  # inputs:
  # linda_data - result of linda()$output
  # intercept - strings corresponding to the name of the intercept
  # outputs:
  # linda_data - renamed data
  
  if (is.data.frame(linda_data)){
    # renaming intercept
    colnames(linda_data) <- gsub("Intercept",group[1],colnames(linda_data))

    colnames(linda_data)[!(grepl(group[1],colnames(linda_data)))] <- paste(group[1],"vs",colnames(linda_data)[!(grepl(group[1],colnames(linda_data)))])
  } else{
    
    # renaming main effects
    names(linda_data)[grepl("Group",names(linda_data)) & !(grepl(":",names(linda_data)))] <-  paste(group[1],"vs",names(linda_data)[grepl("Group",names(linda_data)) & !(grepl(":",names(linda_data)))])
    names(linda_data)[grepl("Country",names(linda_data)) & !(grepl(":",names(linda_data)))] <- paste(group[1], ",", group[2], "-", "CZ vs NO")
    
    # renaming interactions
    names(linda_data)[grepl(":",names(linda_data))] <- paste(group[1],"vs",names(linda_data)[grepl(":",names(linda_data))])
  }
  return(linda_data)
}

binomial_statistics <- function(uni_data, uni_metadata, raw_linda_results,group=NULL,segment){
  # CLR statistics
  data_clr <- vegan::decostand(uni_data,method = "clr", MARGIN = 2, pseudocount=0.5) %>% as.matrix()
  #data_clr <- as.data.frame(apply(uni_data, 2, function(x) x / sum(x)))
  groups <- unique(uni_metadata$Group)
  if (!is.null(group)) groups <- group
  # group1
  group1 <- data_clr[,colnames(data_clr) %in% rownames(uni_metadata)[uni_metadata$Group == groups[1]]]
  res <- apply(group1, 1, function(x) {
    quantile(x, probs = c(0.25,0.5,0.75))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1",groups[1]),paste("MEDIAN",groups[1]),paste("Q3",groups[1]))) %>% 
    as.data.frame()
  
  group2 <- data_clr[,colnames(data_clr) %in% rownames(uni_metadata)[uni_metadata$Group == groups[2]]]
  res <- cbind(res,apply(group2, 1, function(x) {
    quantile(x, probs = c(0.25,0.5,0.75))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1",groups[2]),paste("MEDIAN",groups[2]),paste("Q3",groups[2]))) %>% 
    as.data.frame())
  
  
  # reordering
  #res <- res[,c(grep("Q1",colnames(res)),grep("MEDIAN",colnames(res)),grep("Q3",colnames(res)))]
  if (class(raw_linda_results[[segment]])=="data.frame"){
    raw_linda_result <- raw_linda_results[[segment]]
  } else{
    raw_linda_result <- (raw_linda_results[[segment]][[paste0(groups[1]," vs ","Group",groups[2])]])
  }
  if (is.null(raw_linda_result$Taxonomy)) t <- NA
  else t <- raw_linda_result$Taxonomy
  res <- cbind(ASV=raw_linda_result$ASV,
               Taxonomy=t,
               log2FoldChange=raw_linda_result$log2FoldChange,
               padj=raw_linda_result$padj,
               res)
  res["MEDIAN ALL"] <- apply(data_clr,1,median)
  res["Cliffs_delta"] <- cliffs_delta(data_clr,uni_metadata, group) 
  
  
  # save the results
  if (class(raw_linda_results[[segment]])=="data.frame"){
    raw_linda_results[[segment]] <- res
  } else {
    raw_linda_results[[segment]][[paste0(groups[1]," vs ","Group",groups[2])]] <- res
  }
  return(raw_linda_results)
}

cliffs_delta <- function(data,metadata, group){
  groups <- group
  group1_data <- data[,metadata$Group == groups[1]]
  group2_data <- data[,metadata$Group == groups[2]]
  print(groups[1])
  print(groups[2])
  
  cliffs_deltas <- c()
  for (i in 1:nrow(group1_data)){
    comparison_matrix <- outer(as.matrix(group1_data[i,]), as.matrix(group2_data[i,]), 
                               function(x, y) {ifelse(y > x, 1, -1)}) %>% as.data.frame()
    cd <- sum(comparison_matrix)/(dim(comparison_matrix)[1]*dim(comparison_matrix)[2])
    cliffs_deltas <- c(cliffs_deltas,cd)
  }
  return(cliffs_deltas)
  
}


basic_univariate_statistics <- function(uni_data, group=NULL){
  asv_table <- uni_data[[1]]
  taxa_table <- uni_data[[2]]
  metadata <- uni_data[[3]]
  
  # CLR statistics
  data_clr <- vegan::decostand(asv_table,method = "clr", MARGIN = 2, pseudocount=0.5) %>% as.matrix()
  groups <- unique(metadata$Group)
  if (!is.null(group)) groups <- group
  # group1
  group1 <- data_clr[,colnames(data_clr) %in% rownames(metadata)[metadata$Group == groups[1]]]
  res <- apply(group1, 1, function(x) {
    c(quantile(x, probs = c(0.25, 0.5, 0.75)), mean = mean(x, na.rm = TRUE))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1_clr",groups[1]),
                   paste("MEDIAN_clr",groups[1]),
                   paste("Q3_clr",groups[1]),
                   paste("MEAN_clr",groups[1]))) %>% 
    as.data.frame()
  
  group2 <- data_clr[,colnames(data_clr) %in% rownames(metadata)[metadata$Group == groups[2]]]
  res <- cbind(res,apply(group2, 1, function(x) {
    c(quantile(x, probs = c(0.25, 0.5, 0.75)), 
      mean = mean(x, na.rm = TRUE))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1_clr",groups[2]),
                   paste("MEDIAN_clr",groups[2]),
                   paste("Q3_clr",groups[2]),
                   paste("MEAN_clr",groups[2]))) %>% 
    as.data.frame())
  
  res["MEDIAN_clr_ALL"] <- apply(data_clr,1,median)
  res["Cliffs_delta_clr"] <- cliffs_delta(data_clr,metadata, group) 
  res["LFC_clr"] <- res[,paste("MEAN_clr",groups[2])] - res[,paste("MEAN_clr",groups[1])]

  # RELATIVE ABUNDANCES
  data_ra <- apply(asv_table, 2, function(x) x / sum(x)) %>% as.matrix()
  
  # group1
  group1 <- data_ra[,colnames(data_ra) %in% rownames(metadata)[metadata$Group == groups[1]]]
  res <- cbind(res,apply(group1, 1, function(x) {
    c(quantile(x, probs = c(0.25, 0.5, 0.75)), mean = mean(x, na.rm = TRUE))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1_ra",groups[1]),
                   paste("MEDIAN_ra",groups[1]),
                   paste("Q3_ra",groups[1]),
                   paste("MEAN_ra",groups[1]))) %>% 
    as.data.frame())
  
  group2 <- data_ra[,colnames(data_ra) %in% rownames(metadata)[metadata$Group == groups[2]]]
  res <- cbind(res,apply(group2, 1, function(x) {
    c(quantile(x, probs = c(0.25, 0.5, 0.75)), mean = mean(x, na.rm = TRUE))
  }) %>% t() %>% 
    `colnames<-`(c(paste("Q1_ra",groups[2]),
                   paste("MEDIAN_ra",groups[2]),
                   paste("Q3_ra",groups[2]),
                   paste("MEAN_ra",groups[2]))) %>% 
    as.data.frame())
  
  res["MEDIAN_ra_ALL"] <- apply(data_ra,1,median)
  res["Cliffs_delta_ra"] <- cliffs_delta(data_ra,metadata, group) 
  res["LFC_ra"] <- log((res[,paste("MEAN_ra",groups[2])]) / 
                         (res[,paste("MEAN_ra",groups[1])]),2)
  res[paste("PREVALENCE_",groups[1])] <- (rowSums(group1>0))/(ncol(group1))
  res[paste("PREVALENCE_",groups[2])] <- (rowSums(group2>0))/(ncol(group2))
  
  # reordering
  res <- res[,c(grep("Q1",colnames(res)),
                grep("MEDIAN",colnames(res)),
                grep("Q3",colnames(res)),
                grep("MEAN",colnames(res)),
                grep("Cliffs_delta",colnames(res)),
                grep("LFC",colnames(res)),
                grep("PREVALENCE",colnames(res)))]
  
  # taxonomy
  res %<>% rownames_to_column("SeqID")
  #res <- merge(res %>% rownames_to_column("SeqID"),taxa_table,by="SeqID",all = TRUE)
  
  return(res)
}


## glmnet functions ----
glmnet_binomial <- function(data,
                            outcome,
                            sample_method = 'atypboot',
                            clust_var=NULL,
                            N = 10, # number of bootstrap datasets
                            alphas = seq(0, 1, by = 0.2),
                            family = 'binomial',
                            overfitting_check=FALSE,
                            seed = 123,
                            reuse=FALSE,
                            file=NULL,
                            Q=NULL) {
  
  
  logit <- function(x){log(x/(1-x))}
  inv_logit <- function(x){exp(x)/(1+exp(x))}
  if (all(data[,1] >= 0 & data[,1] <= 1)) ra = TRUE
  else ra = FALSE
  
  if (overfitting_check) {
    data$Group <- sample(data$Group)
  }
  
  if (reuse){
    if (ra) {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"enet_model_ra.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"enet_model_ra.RData"))
    } else {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"enet_model.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"enet_model.RData"))
    }
  } else {
    set.seed(seed)
    
    ## where to save relevant information
    auc_validation <- vector('double', N)
    accuracy_validation <- vector('double', N)
    predictions <- vector("list", N)
    conf_matrices <- vector("list", N)
    kfold_rocobjs <- vector("list", N)
    
    ## original data in a matrix form
    original_outcome <- base::as.matrix(data[[outcome]])
    
    if(!is.numeric(original_outcome)){
      original_outcome <- factor(original_outcome)
      mapping <- setNames(levels(original_outcome), 0:(length(levels(original_outcome))-1))
      original_outcome <- as.numeric(factor(original_outcome)) - 1 
    }
    
    original_predictors <- data %>% 
      dplyr::select(-dplyr::all_of(c(outcome,"Country",clust_var))) %>% 
      as.matrix()
    
    # prediction on original sample
    
    ## optimize lambda and alpha
    lamb_1se <- vector('double', length(alphas))
    alpha <- vector('double', length(alphas))
    deviance <- vector('double', length(alphas))
    
    for(a in seq_along(alphas)){
      tr <- cv.glmnet(x = original_predictors, 
                      y = original_outcome, 
                      alpha = alphas[a], 
                      family = family,
                      type.measure = 'deviance')
      lamb_1se[a] <- tr[["lambda.1se"]]
      
      alpha[a] = alphas[a]
      deviance[a] = tr$cvm[which(tr$lambda == tr[["lambda.1se"]])]
    }
    
    optim_par <- data.frame(lamb_1se, alpha, deviance) %>% 
      arrange(deviance)
    
    
    ## fit with optimized hyperparameters
    fit <- glmnet(x = original_predictors, 
                  y = original_outcome, 
                  alpha = optim_par$alpha[1],
                  lambda = optim_par$lamb_1se[1],
                  family = family)
    
    
    predicted_orig = as.numeric(predict(fit, newx = original_predictors,type = "class"))
    predicted_num = as.numeric(predict(fit, newx = original_predictors))
    
    
    ## get predictions and performance
    prediction <- data.frame(
      predicted = predicted_orig,
      predicted_num=predicted_num,
      outcome = original_outcome)
    
    if (length(unique(prediction$predicted))>1){
      conf_matrix_orig <- table(True = prediction$outcome, Predicted = prediction$predicted) 
    } else {
      conf_matrix_orig <- cbind(table(True = prediction$outcome, Predicted = prediction$predicted),c(0,0)) 
    }
    
    rocobj <- roc(outcome ~ predicted_num, 
                  data = prediction,
                  direction = '<',
                  levels = c(0, 1))
    
    # country specific  - czech
    if ("CZ" %in% data$Country){ 
      ## predicted classes as numerical variables
      predictors_czech <- data[data$Country == "CZ",] %>% 
        dplyr::select(-dplyr::all_of(c(outcome,"Country", clust_var))) %>% 
        as.matrix()
      
      predicted_czech = as.numeric(predict(fit, newx = predictors_czech,type = "class"))
      predicted_czech_num = as.numeric(predict(fit, newx = predictors_czech))
      outcome_czech <- original_outcome[data$Country == "CZ"]
      
      ## get predictions and performance
      prediction_czech <- data.frame(
        predicted = predicted_czech,
        predicted_num=predicted_czech_num,
        outcome = outcome_czech)
      
      ## confusion matrix
      if (length(unique(prediction_czech$predicted))>1){
        conf_matrix_czech <- table(True = prediction_czech$outcome, Predicted = prediction_czech$predicted) 
      } else {
        conf_matrix_czech <- cbind(table(True = prediction_czech$outcome, Predicted = prediction_czech$predicted),c(0,0)) 
      }
      
      auc_czech <- roc(outcome ~ predicted_num, 
                       data = prediction_czech,
                       direction = '<',
                       levels = c(0, 1))$auc
    } else{
      prediction_czech = NULL
      auc_czech <- NaN
      conf_matrix_czech = NaN
    }
    # country specific  - norway
    ## predicted classes as numerical variables
    if ("NO" %in% data$Country){
      predictors_no <- data[data$Country == "NO",] %>% 
        dplyr::select(-dplyr::all_of(c(outcome,"Country",clust_var))) %>% 
        as.matrix()
      
      predicted_no = as.numeric(predict(fit, newx = predictors_no,type = "class"))
      predicted_no_num = as.numeric(predict(fit, newx = predictors_no))
      outcome_no <- original_outcome[data$Country == "NO"]
      
      ## get predictions and performance
      prediction_no <- data.frame(
        predicted = predicted_no,
        predicted_num=predicted_no_num,
        outcome = outcome_no)
      
      ## confusion matrix
      if (length(unique(prediction_no$predicted))>1){
        conf_matrix_no <- table(True = prediction_no$outcome, Predicted = prediction_no$predicted) 
      } else {
        conf_matrix_no <- cbind(table(True = prediction_no$outcome, Predicted = prediction_no$predicted),c(0,0)) 
      }
      
      auc_no <- roc(outcome ~ predicted_num, 
                    data = prediction_no,
                    direction = '<',
                    levels = c(0, 1))$auc
    } else{
      prediction_no = NULL
      auc_no <- NaN
      conf_matrix_no = NaN
    }
    
    fitted <- data.frame(
      alpha = optim_par$alpha[1],
      lambda = optim_par$lamb_1se[1],
      auc = rocobj$auc,
      auc_czech=auc_czech,
      auc_no=auc_no,
      accuracy = mean(prediction$predicted== prediction$outcome),
      accuracy_num=mean(ifelse(
        prediction$predicted > 0, 1, 0) == prediction$outcome),
      accuracy_czech=mean(prediction_czech$preebnndicted== prediction_czech$outcome),
      accuracy_no=mean(prediction_no$predicted== prediction_no$outcome))
    
    
    ## simulated data
  
    sampled_data <- sampler(data, 
                            outcome = outcome,
                            clust_var = clust_var,
                            sample_method = sample_method,
                            N = N, 
                            seed = seed)
    
    if(!is.null(clust_var)) clust_var <- "id"
    
    for (i in 1:N){
      
      ## sampled data in a matrix form
      sampled_outcome <- as.matrix(sampled_data[[1]][[i]]$outcome)
      if(!is.numeric(sampled_outcome)){
        sampled_outcome <- factor(sampled_outcome)
        mapping <- setNames(levels(sampled_outcome), 0:(length(levels(sampled_outcome))-1))
        sampled_outcome <- as.numeric(factor(sampled_outcome)) - 1 
      }
      
      sampled_predictors <- sampled_data[[1]][[i]] %>% 
        dplyr::select(-dplyr::all_of(c("outcome","id", "obs_id","Country", clust_var))) %>%
        as.matrix()
      
      ## re-optimize alpha and lambda
      lamb_1se <- vector('double', length(alphas))
      alpha <- vector('double', length(alphas))
      deviance <- vector('double', length(alphas))
      
      for(a in seq_along(alphas)){
        tr <- cv.glmnet(x = sampled_predictors, 
                        y = sampled_outcome, 
                        alpha = alphas[a], 
                        family = family,
                        type.measure = 'deviance')
        lamb_1se[a] <- tr[["lambda.1se"]]
        
        alpha[a] = alphas[a]
        deviance[a] = tr$cvm[which(tr$lambda == tr[["lambda.1se"]])]
      }
      
      optim_par <- data.frame(lamb_1se, alpha, deviance) %>% 
        arrange(deviance)
      
      
      ## fit models with re-optimized hyperparameters
      sampled_fit <- glmnet(sampled_predictors,
                            sampled_outcome, 
                            alpha = optim_par$alpha[1],
                            lambda = optim_par$lamb_1se[1],
                            family = family)
      
      
      valid_outcome <- as.matrix(sampled_data[[2]][[i]]$outcome)
      valid_predictors <- sampled_data[[2]][[i]] %>% 
        dplyr::select(-outcome, -id, - Country, -obs_id) %>% 
        as.matrix()
      
      prediction_onValidation <- data.frame(
        predicted = as.numeric(predict(sampled_fit, newx = valid_predictors,type = "class")),
        predicted_num=as.numeric(predict(sampled_fit, newx = valid_predictors)),
        outcome = valid_outcome,
        iteration = i)
      
      ## record performances
      kfold_rocobjs[[i]] <- roc(outcome ~ predicted_num, 
                              data = prediction_onValidation,
                              direction = '<',
                              levels = c(0, 1))
      
      auc_validation[i] <- kfold_rocobjs[[i]]$auc
      
      accuracy_validation[i] <-  mean(prediction_onValidation$predicted == prediction_onValidation$outcome)
      predictions[[i]] <- data.frame(prediction_onValidation, id = sampled_data[[2]][[i]]$id)
      
      if (length(unique(prediction_onValidation$predicted))>1){
        conf_matrices[[i]] <- table(True = prediction_onValidation$outcome, Predicted = prediction_onValidation$predicted) 
      } else {
        conf_matrices[[i]] <- cbind(table(True = prediction_onValidation$outcome, Predicted = prediction_onValidation$predicted),c(0,0)) 
      }
    }
    
    ## connect predictions
    predictions <- bind_rows(predictions)
    
    predictions2 <- predictions %>% 
      group_by(id) %>%
      summarise(
        predicted = mean(predicted_num),
        outcome = mean(outcome)
      ) %>% 
      ungroup()
    
    conf_matrices_list <- list(original=conf_matrix_orig,
                               czech=conf_matrix_czech,
                               no=conf_matrix_no)
    
    ## aggregate obtained information
    valid_performances <- data.frame(
      auc_validation,
      accuracy_validation)
    
    betas <- fit$beta
    
    model_summary <- fitted %>% 
      mutate(
        auc_optimism_corrected = mean(valid_performances$auc_validation),
        auc_optimism_corrected_CIL = quantile(valid_performances$auc_validation, probs = 0.025),
        auc_optimism_corrected_CIU = quantile(valid_performances$auc_validation, probs = 0.975),
        accuracy_optimism_corrected = mean(valid_performances$accuracy_validation),
        accuracy_optimism_corrected_CIL =  quantile(valid_performances$accuracy_validation, probs = 0.025),
        accuracy_optimism_corrected_CIU =  quantile(valid_performances$accuracy_validation, probs = 0.975)) %>% 
      dplyr::select(alpha, lambda, 
                    auc, auc_czech, auc_no,  
                    auc_optimism_corrected:auc_optimism_corrected_CIU,
                    accuracy, accuracy_czech, accuracy_no, 
                    accuracy_optimism_corrected:accuracy_optimism_corrected_CIU)
    
    # calibration curve
    calibration_plot <- suppressWarnings(
      predictions2 %>% 
        mutate(iteration = factor('A'),
               predicted = inv_logit(predicted)) %>% 
        ggplot(aes(x = predicted, y = outcome, group = iteration)) +
        geom_smooth(data = predictions,
                    aes(x = inv_logit(predicted_num), y = outcome),
                    se = FALSE,
                    color = 'grey35',
                    linewidth = 0.1,
                    method = 'loess',
                    span = 2/log10(nrow(prediction))) +
        geom_smooth(method = 'loess',
                    se = TRUE,
                    color = 'red',
                    fill = 'red', 
                    alpha = 0.25, 
                    span = 2/log10(nrow(prediction))) + 
        coord_cartesian(x = c(min(inv_logit(predictions$predicted_num)),
                              max(inv_logit(predictions$predicted_num))), 
                        y = c(0,1)) + 
        geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = 'dashed') +
        labs(x = "Prediction", y = "Outcome") + theme_minimal()
    )
    
    ## define outputs
    enet_model <- list(model_summary = model_summary, 
                       valid_performances = valid_performances, 
                       predictions = prediction, 
                       betas = betas,
                       conf_matrices=conf_matrices_list,
                       rocobj=rocobj,
                       kfold_rocobjs=kfold_rocobjs,
                       plot = calibration_plot, 
                       trained_model=fit)
    
    # save results
    if (overfitting_check){
     if (!dir.exists(file.path("../intermediate_files/models_overfitting_check/",Q,file))){
      dir.create(file.path("../intermediate_files/models_overfitting_check/",Q,file))
     }
      if (ra) save(enet_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"enet_model_ra.RData"))
      else save(enet_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"enet_model.RData"))
      
    } else {
      if (!dir.exists(file.path("../intermediate_files/models/",Q,file))){
      dir.create(file.path("../intermediate_files/models/",Q,file))
      }
      if (ra) save(enet_model,file=file.path("../intermediate_files/models/",Q,file,"enet_model_ra.RData"))
      else save(enet_model,file=file.path("../intermediate_files/models/",Q,file,"enet_model.RData"))
    }
    
  }
  
  return(enet_model)
}


sampler <- function(data, outcome,
                    seed = NULL, 
                    clust_var=NULL,
                    sample_method = 'atypboot',
                    N = 10){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  data <- data %>% 
    mutate(obs_id = as.character(1:nrow(data)))
  
  if (is.null(clust_var)) clust_var="obs_id" 
  if(clust_var != 'id'){
    data <- data %>% 
      mutate(id = data[[clust_var]]) %>% 
      dplyr::select(-dplyr::all_of(clust_var))
  }
  
  if(colnames(data[outcome]) != 'outcome'){
    data <- data %>% 
      mutate(outcome = data[[outcome]]) %>% 
      dplyr::select(-dplyr::all_of(c(outcome)))
  }
  
  data <- data %>% 
    mutate(obs_id = as.character(1:nrow(data))) %>%
    dplyr::select(id,dplyr::everything())
  
 if (sample_method == 'atypboot') {
    unique_ids <- unique(data$id)
    id_groups <- sample(rep(1:N, length.out = length(unique_ids)))
    
    train_data <- list()
    valid_data <- list()
    
    for (i in 1:N) {
      repeat {
        train <- data.frame(id = sample(unique(data$id), 
                                        length(unique(data$id)), 
                                        replace = TRUE))
        
        temp_train <- train %>% 
          left_join(data, by = 'id', relationship = "many-to-many")
        
        temp_valid <- data.frame(
          obs_id = data[!data$obs_id %in% temp_train$obs_id, ]$obs_id) %>%
          left_join(data, by = 'obs_id')
        
        if (!mean(temp_train$outcome) %in% c(0, 1) & !mean(temp_valid$outcome) %in% c(0, 1)) {
          train_data[[i]] <- temp_train
          valid_data[[i]] <- temp_valid
          {break} 
        }
      }
    }
    return(list(train_data, valid_data))
  }

}

## Machine learning functions ----
gbm_binomial <- function(data,
                         outcome,
                         sample_method = 'atypboot',
                         clust_var=NULL,
                         N = 10, # number of bootstrap datasets
                         family = 'binomial',
                         overfitting_check=FALSE,
                         seed = 123,
                         reuse=FALSE,
                         file=NULL,
                         Q=NULL) {
  
  #n.trees, interaction.depth, shrinkage, n.minobsinnode
  if (all(data[,1] >= 0 & data[,1] <= 1)) ra = TRUE
  else  ra = FALSE
  
  if (overfitting_check) {
    data$Group <- sample(data$Group)
  }
  
  if (reuse){
    if (ra) {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"gbm_model_ra.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"gbm_model_ra.RData"))
    } else {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"gbm_model.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"gbm_model.RData"))
    }
  } else {
    set.seed(seed)
    
    data <- data.frame(data,check.names=TRUE)
    
    ## where to save relevant information
    auc_validation <- vector('double', N)
    accuracy_validation <- vector('double', N)
    predictions <- vector("list", N)
    conf_matrices <- vector("list", N)
    kfold_rocobjs <- vector("list", N)
    
    ## original data in a matrix form
    original_outcome <- base::as.matrix(data[[outcome]])
    
    if(!is.numeric(original_outcome)){
      original_outcome <- factor(original_outcome)
      mapping <- setNames(levels(original_outcome), 0:(length(levels(original_outcome))-1))
      original_outcome <- as.numeric(factor(original_outcome)) - 1 
    }
    
    original_predictors <- data %>% 
      dplyr::select(-dplyr::all_of(c(outcome,"Country",clust_var))) %>% 
      as.matrix()
    
    original_outcome <- as.factor(original_outcome)
    
    tune_grid <- expand.grid(
      n.trees = c(100,200,500),  # Number of predictors to consider
      interaction.depth = c(1,3,5),         # Minimum node size
      shrinkage = 0.1,
      n.minobsinnode = c(10,20,30)
    )
    
    fitControl <- trainControl(method = "repeatedcv",
                               ## 5-fold CV...
                               number = 5,
                               ## repeated 5 times
                               repeats = 1)
    
    # OPTIMALIZATION Train the random forest model using ranger
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
    
    gbm_model <- caret::train(
      x = original_predictors, 
      y = original_outcome, 
      method = "gbm",  # Use kknn method for more flexibility
      tuneGrid = tune_grid,  # Manhattan or Euclidean distance
      metric = "Accuracy", # Specify the evaluation metric
      trControl = fitControl
    )
    
    stopCluster(cl)
    
    optim_n.trees = gbm_model$bestTune$n.trees
    optim_interaction.depth = gbm_model$bestTune$interaction.depth
    optim_shrinkage = gbm_model$bestTune$shrinkage
    optim_n.minobsinnode = gbm_model$bestTune$n.minobsinnode
    
    # Combine predictors and outcome into a single data frame
    original_outcome <- as.numeric(original_outcome) - 1
    train_data <- data.frame(original_outcome, 
                             original_predictors, 
                             check.names = TRUE)
    
    # Correct formula: specify y explicitly
    fit <- gbm(
      formula = original_outcome ~ .,         
      distribution = "bernoulli",
      data = train_data,        # Use the combined data frame
      n.trees = optim_n.trees,  # Use the best value of k from tuning
      interaction.depth = optim_interaction.depth,
      shrinkage = optim_shrinkage,
      n.minobsinnode = optim_n.minobsinnode
    )
    
    predicted_num = predict(fit, newdata = as.data.frame(original_predictors),type="response")
    predicted_orig = ifelse(predicted_num > 0.5, 1, 0)
    
    prediction <- data.frame(
      predicted = predicted_orig,
      predicted_num=predicted_num,
      outcome = original_outcome)
    
    rocobj <- roc(outcome ~ predicted_num, 
                  data = prediction,
                  direction = '<',
                  levels = c(0, 1))
    
    fitted <- data.frame(
      n.trees = optim_n.trees,
      interaction.depth = optim_interaction.depth,
      shrinkage = optim_shrinkage,
      n.minobsinnode = optim_n.minobsinnode,
      auc = rocobj$auc,
      accuracy = mean(prediction$predicted==prediction$outcome),
      accuracy_num=mean(ifelse(prediction$predicted_num > 0.5, 1, 0) == prediction$outcome)
    )
    
    # BOOTSTRAP
    sampled_data <- sampler(data, 
                            outcome = outcome,
                            clust_var = clust_var,
                            sample_method = sample_method,
                            N = N, 
                            seed = seed)
    
    if(!is.null(clust_var)) clust_var <- "id"
    
    for (i in 1:N){
      
      ## sampled data in a matrix form
      sampled_outcome <- as.matrix(sampled_data[[1]][[i]]$outcome)
      
      if (!is.numeric(sampled_outcome)){
        sampled_outcome <- factor(sampled_outcome)
        mapping <- setNames(levels(sampled_outcome), 0:(length(levels(sampled_outcome))-1))
        sampled_outcome <- as.numeric(factor(sampled_outcome)) - 1 
      }
      
      sampled_outcome <- as.factor(sampled_outcome)
      
      sampled_predictors <- sampled_data[[1]][[i]] %>% 
        dplyr::select(-dplyr::all_of(c("outcome","id", "obs_id","Country", clust_var))) %>%
        as.matrix()
      
      ## re-optimize parameters
      # OPTIMALIZATION
      
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      
      gbm_model_sampled <- caret::train(
        x = sampled_predictors, 
        y = sampled_outcome, 
        method = "gbm",  # Use kknn method for more flexibility
        tuneGrid = tune_grid,  # Manhattan or Euclidean distance
        metric = "Accuracy", # Specify the evaluation metric
        trControl = fitControl
      )
      
      stopCluster(cl)
      
      
      optim_n.trees_sampled = gbm_model_sampled$bestTune$n.trees
      optim_interaction.depth_sampled = gbm_model_sampled$bestTune$interaction.depth
      optim_shrinkage_sampled = gbm_model_sampled$bestTune$shrinkage
      optim_n.minobsinnode_sampled = gbm_model_sampled$bestTune$n.minobsinnode
      
      # Combine predictors and outcome into a single data frame
      sampled_outcome <- as.numeric(sampled_outcome) - 1
      train_data_sampled <- data.frame(sampled_outcome, 
                                       sampled_predictors, 
                                       check.names = TRUE)
      # TRAINING
      sampled_fit <- gbm(
        formula = sampled_outcome ~ .,         
        distribution = "bernoulli",
        data = train_data_sampled,        # Use the combined data frame
        n.trees = optim_n.trees,  # Use the best value of k from tuning
        interaction.depth = optim_interaction.depth,
        shrinkage = optim_shrinkage,
        n.minobsinnode = optim_n.minobsinnode
      )
      
      sampled_predicted_num = predict(sampled_fit, newdata = as.data.frame(sampled_predictors),type="response")
      sampled_predicted_orig = ifelse(sampled_predicted_num > 0.5, 1, 0)
      
      
      # VALIDATION
      valid_outcome <- as.factor(as.matrix(sampled_data[[2]][[i]]$outcome))
      valid_predictors <- sampled_data[[2]][[i]] %>% 
        dplyr::select(-outcome, -id, - Country, -obs_id) %>% 
        as.matrix()
      
      valid_predicted_num = predict(sampled_fit, 
                                    newdata = as.data.frame(valid_predictors))
      
      valid_predicted_orig <- ifelse(valid_predicted_num > 0.5, 1, 0)
      
      
      prediction_onValidation <- data.frame(
        predicted = valid_predicted_orig,
        predicted_num = valid_predicted_num,
        outcome = valid_outcome,
        iteration = i)
      
      ## record performances
      kfold_rocobjs[[i]] <- roc(outcome ~ predicted_num, 
                                data = prediction_onValidation,
                                direction = '<',
                                levels = c(0, 1))
      
      auc_validation[i] <- kfold_rocobjs[[i]]$auc
      accuracy_validation[i] <- mean(prediction_onValidation$predicted == prediction_onValidation$outcome)
      predictions[[i]] <- data.frame(prediction_onValidation, id = sampled_data[[2]][[i]]$id)
    }
    
    ## connect predictions
    predictions <- bind_rows(predictions)
    
    ## aggregate obtained information
    valid_performances <- data.frame(
      auc_validation,
      accuracy_validation)
    
    model_summary <- fitted %>% 
      mutate(
        auc_optimism_corrected = mean(valid_performances$auc_validation),
        auc_optimism_corrected_CIL = quantile(valid_performances$auc_validation, probs = 0.025),
        auc_optimism_corrected_CIU = quantile(valid_performances$auc_validation, probs = 0.975),
        accuracy_optimism_corrected = mean(valid_performances$accuracy_validation),
        accuracy_optimism_corrected_CIL =  quantile(valid_performances$accuracy_validation, probs = 0.025),
        accuracy_optimism_corrected_CIU =  quantile(valid_performances$accuracy_validation, probs = 0.975)) %>% 
      dplyr::select(n.trees,
                    interaction.depth,
                    shrinkage,
                    n.minobsinnode,
                    auc,  
                    auc_optimism_corrected:auc_optimism_corrected_CIU,
                    accuracy,
                    accuracy_optimism_corrected:accuracy_optimism_corrected_CIU)
    
    ## define outputs
    gbm_model <- list(model_summary = model_summary, 
                     valid_performances = valid_performances, 
                     predictions = prediction, 
                     rocobj=rocobj,
                     kfold_rocobjs=kfold_rocobjs,
                     trained_model=fit)
    
    # save results
    if (overfitting_check){
      if (!dir.exists(file.path("../intermediate_files/models_overfitting_check/",Q,file))){
        dir.create(file.path("../intermediate_files/models_overfitting_check/",Q,file))
      } 
      if (ra) save(gbm_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"gbm_model_ra.RData"))
      else save(gbm_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"gbm_model.RData"))
    } else {
      if (!dir.exists(file.path("../intermediate_files/models/",Q,file))){
        dir.create(file.path("../intermediate_files/models/",Q,file))
      }
      if (ra) save(gbm_model,file=file.path("../intermediate_files/models/",Q,file,"gbm_model_ra.RData"))
      else save(gbm_model,file=file.path("../intermediate_files/models/",Q,file,"gbm_model.RData"))
    }
  }
  return(gbm_model)
}

knn_binomial <- function(data,
                         outcome,
                         sample_method = 'atypboot',
                         clust_var=NULL,
                         N = 10, # number of bootstrap datasets
                         family = 'binomial',
                         overfitting_check=FALSE,
                         seed = 123,
                         reuse=FALSE,
                         file=NULL,
                         Q=NULL) {
  
  if (all(data[,1] >= 0 & data[,1] <= 1)) ra = TRUE
  else  ra = FALSE
  
  if (overfitting_check) {
    data$Group <- sample(data$Group)
  }
  
  if (reuse){
    if (ra) {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"knn_model_ra.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"knn_model_ra.RData"))
    } else {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"knn_model.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"knn_model.RData"))
    }
  } else {
    set.seed(seed)
    
    data <- data.frame(data,check.names=TRUE)
    
    ## where to save relevant information
    auc_validation <- vector('double', N)
    accuracy_validation <- vector('double', N)
    predictions <- vector("list", N)
    conf_matrices <- vector("list", N)
    kfold_rocobjs <- vector("list", N)
    
    ## original data in a matrix form
    original_outcome <- base::as.matrix(data[[outcome]])
    
    if(!is.numeric(original_outcome)){
      original_outcome <- factor(original_outcome)
      mapping <- setNames(levels(original_outcome), 0:(length(levels(original_outcome))-1))
      original_outcome <- as.numeric(factor(original_outcome)) - 1 
    }
    
    original_predictors <- data %>% 
      dplyr::select(-dplyr::all_of(c(outcome,"Country",clust_var))) %>% 
      as.matrix()
    
    tune_grid <- expand.grid(k = seq(10, 30, by = 1))
    
    fitControl <- trainControl(method = "repeatedcv",
                               ## 5-fold CV...
                               number = 5,
                               ## repeated 5 times
                               repeats = 5)
    
    # OPTIMALIZATION
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)
    
    knn_model <- caret::train(
      x = original_predictors, 
      y = as.factor(original_outcome), 
      method = "knn",  # Use kknn method for more flexibility
      tuneGrid = tune_grid,  # Manhattan or Euclidean distance
      metric = "Accuracy", # Specify the evaluation metric
      trControl = fitControl
    )
    
    stopCluster(cl)
    optim_k = knn_model$bestTune$k
    
    # Combine predictors and outcome into a single data frame
    train_data <- data.frame(original_outcome, 
                             original_predictors, 
                             check.names = TRUE)
    
    # Correct formula: specify y explicitly
    fit <- knn3(
      formula = original_outcome ~ .,          # Explicitly specify the outcome variable 'y'
      data = train_data,        # Use the combined data frame
      k = optim_k  # Use the best value of k from tuning
    )
    
    predicted_num = predict(fit, 
                            newdata = as.data.frame(original_predictors))[,2]
    
    predicted_orig <- ifelse(predicted_num > 0.5, 1, 0)
    
    prediction <- data.frame(
      predicted = predicted_orig,
      predicted_num=predicted_num,
      outcome = original_outcome)
    
    rocobj <- roc(outcome ~ predicted_num, 
                  data = prediction,
                  direction = '<',
                  levels = c(0, 1))
    
    fitted <- data.frame(
      k = optim_k,
      auc = rocobj$auc,
      accuracy = mean(prediction$predicted==prediction$outcome),
      accuracy_num=mean(ifelse(prediction$predicted_num > 0.5, 1, 0) == prediction$outcome))
  
    # BOOTSTRAP
    sampled_data <- sampler(data, 
                            outcome = outcome,
                            clust_var = clust_var,
                            sample_method = sample_method,
                            N = N, 
                            seed = seed)
    
    if(!is.null(clust_var)) clust_var <- "id"
    
    for (i in 1:N){
      
      ## sampled data in a matrix form
      sampled_outcome <- as.matrix(sampled_data[[1]][[i]]$outcome)
      
      if(!is.numeric(sampled_outcome)){
        sampled_outcome <- factor(sampled_outcome)
        mapping <- setNames(levels(sampled_outcome), 0:(length(levels(sampled_outcome))-1))
        sampled_outcome <- as.numeric(factor(sampled_outcome)) - 1 
      }
      
      sampled_predictors <- sampled_data[[1]][[i]] %>% 
        dplyr::select(-dplyr::all_of(c("outcome","id", "obs_id","Country", clust_var))) %>%
        as.matrix()
      
      ## re-optimize k
      # OPTIMALIZATION
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      
      knn_model_sampled <- caret::train(
        x = sampled_predictors, 
        y = as.factor(sampled_outcome), 
        method = "knn",  # Use kknn method for more flexibility
        tuneGrid = tune_grid,  # Manhattan or Euclidean distance
        metric = "Accuracy", # Specify the evaluation metric
        trControl = fitControl
      )
      stopCluster(cl)
      
      optim_k_sampled = knn_model_sampled$bestTune$k
      
      # Combine predictors and outcome into a single data frame
      train_data_sampled <- data.frame(sampled_outcome, 
                               sampled_predictors, 
                               check.names = TRUE)
      
      
      # TRAINING
      sampled_fit <- knn3(
        formula = sampled_outcome ~ .,          # Explicitly specify the outcome variable 'y'
        data = train_data_sampled,        # Use the combined data frame
        k = optim_k_sampled  # Use the best value of k from tuning
      )
      
      # VALIDATION
      valid_outcome <- as.matrix(sampled_data[[2]][[i]]$outcome)
      valid_predictors <- sampled_data[[2]][[i]] %>% 
        dplyr::select(-outcome, -id, - Country, -obs_id) %>% 
        as.matrix()
      
      valid_predicted_num = predict(sampled_fit, 
                              newdata = as.data.frame(valid_predictors),type="prob")[,2]
  
      valid_predicted_orig <- ifelse(valid_predicted_num > 0.5, 1, 0)
      
      prediction_onValidation <- data.frame(
        predicted = valid_predicted_orig,
        predicted_num=valid_predicted_num,
        outcome = valid_outcome,
        iteration = i)
      
      ## record performances
      kfold_rocobjs[[i]] <- roc(outcome ~ predicted_num, 
                                data = prediction_onValidation,
                                direction = '<',
                                levels = c(0, 1))
      
      auc_validation[i] <- kfold_rocobjs[[i]]$auc
      accuracy_validation[i] <- mean(prediction_onValidation$predicted == prediction_onValidation$outcome)
      predictions[[i]] <- data.frame(prediction_onValidation, id = sampled_data[[2]][[i]]$id)
    }
    
    ## connect predictions
    predictions <- bind_rows(predictions)
    
    predictions2 <- predictions %>% 
      group_by(id) %>%
      summarise(
        predicted = mean(predicted_num),
        outcome = mean(outcome)
      ) %>% 
      ungroup()
    
    ## aggregate obtained information
    valid_performances <- data.frame(
      auc_validation,
      accuracy_validation)
    
    model_summary <- fitted %>% 
      mutate(
        auc_optimism_corrected = mean(valid_performances$auc_validation),
        auc_optimism_corrected_CIL = quantile(valid_performances$auc_validation, probs = 0.025),
        auc_optimism_corrected_CIU = quantile(valid_performances$auc_validation, probs = 0.975),
        accuracy_optimism_corrected = mean(valid_performances$accuracy_validation),
        accuracy_optimism_corrected_CIL =  quantile(valid_performances$accuracy_validation, probs = 0.025),
        accuracy_optimism_corrected_CIU =  quantile(valid_performances$accuracy_validation, probs = 0.975)) %>% 
      dplyr::select(k, auc,  
                    auc_optimism_corrected:auc_optimism_corrected_CIU,
                    accuracy,
                    accuracy_optimism_corrected:accuracy_optimism_corrected_CIU)
    
    ## define outputs
    knn_model <- list(model_summary = model_summary, 
                       valid_performances = valid_performances, 
                       predictions = prediction, 
                       rocobj=rocobj,
                       kfold_rocobjs=kfold_rocobjs,
                       trained_model=fit)
    
    # save results
    if (overfitting_check){
      if (!dir.exists(file.path("../intermediate_files/models_overfitting_check/",Q,file))){
        dir.create(file.path("../intermediate_files/models_overfitting_check/",Q,file))
      } 
      if (ra) save(knn_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"knn_model_ra.RData"))
      else save(knn_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"knn_model.RData"))
    } else {
      if (!dir.exists(file.path("../intermediate_files/models/",Q,file))){
        dir.create(file.path("../intermediate_files/models/",Q,file))
      }
      if (ra) save(knn_model,file=file.path("../intermediate_files/models/",Q,file,"knn_model_ra.RData"))
      else save(knn_model,file=file.path("../intermediate_files/models/",Q,file,"knn_model.RData"))
    }
    
  }
  
  return(knn_model)
  
}

rf_binomial <- function(data,
                         outcome,
                         sample_method = 'atypboot',
                         clust_var=NULL,
                         N = 10, # number of bootstrap datasets
                         alphas = seq(0, 1, by = 0.2),
                         family = 'binomial',
                         overfitting_check=FALSE,
                         seed = 123,
                         reuse=FALSE,
                         file=NULL,
                         Q=NULL) {

  # https://topepo.github.io/caret/train-models-by-tag.html#random-forest
  if (all(data[,1] >= 0 & data[,1] <= 1)) ra = TRUE
  else  ra = FALSE
  
  if (overfitting_check) {
    data$Group <- sample(data$Group)
  }
  
  if (reuse){
    if (ra) {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model_ra.RData"))
      load(file.path("../intermediate_files/models/",Q,file,"rf_model_ra.RData"))
    } else {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model.RData"))
      load(file.path("../intermediate_files/models/",Q,file,"rf_model.RData"))
    }
  } else {
    set.seed(seed)
    
    data <- data.frame(data,check.names=TRUE)
    
    ## where to save relevant information
    auc_validation <- vector('double', N)
    accuracy_validation <- vector('double', N)
    predictions <- vector("list", N)
    conf_matrices <- vector("list", N)
    kfold_rocobjs <- vector("list", N)
    
    ## original data in a matrix form
    original_outcome <- base::as.matrix(data[[outcome]])
    
    if(!is.numeric(original_outcome)){
      original_outcome <- factor(original_outcome)
      mapping <- setNames(levels(original_outcome), 0:(length(levels(original_outcome))-1))
      original_outcome <- as.numeric(factor(original_outcome)) - 1 
    }
    
    original_predictors <- data %>% 
      dplyr::select(-dplyr::all_of(c(outcome,"Country",clust_var))) %>% 
      as.matrix()
    
    original_outcome <- as.factor(original_outcome)
    
    tune_grid <- expand.grid(
      mtry = seq(1, ncol(original_predictors),2),  # Number of predictors to consider
      min.node.size = c(2, 5),         # Minimum node size
      splitrule = "gini"
    )
    
    fitControl <- trainControl(method = "repeatedcv",
                               ## 5-fold CV...
                               number = 5,
                               ## repeated 5 times
                               repeats = 1)
    
    # OPTIMALIZATION Train the random forest model using ranger
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)

    rf_model <- caret::train(
      x = original_predictors, 
      y = original_outcome, 
      method = "ranger",  # Use kknn method for more flexibility
      tuneGrid = tune_grid,  # Manhattan or Euclidean distance
      metric = "Accuracy", # Specify the evaluation metric
      trControl = fitControl,
      importance = "impurity"
    )
    
    stopCluster(cl)
    
    optim_mtry = rf_model$bestTune$mtry
    optim_splitrule = rf_model$bestTune$splitrule
    optim_min.node.size = rf_model$bestTune$min.node.size
    
    # Combine predictors and outcome into a single data frame
    train_data <- data.frame(original_outcome, 
                             original_predictors, 
                             check.names = TRUE)
    
    # Correct formula: specify y explicitly
    fit <- ranger(
      formula = original_outcome ~ .,          # Explicitly specify the outcome variable 'y'
      data = train_data,        # Use the combined data frame
      mtry = optim_mtry,  # Use the best value of k from tuning
      splitrule = optim_splitrule,
      min.node.size = optim_min.node.size,
      probability = TRUE
    )
    
    predicted_num = predict(fit, 
                            data = as.data.frame(original_predictors)
                            )$predictions[,2]
    
    predicted_orig <- ifelse(predicted_num > 0.5, 1, 0)
    
    prediction <- data.frame(
      predicted = predicted_orig,
      predicted_num=predicted_num,
      outcome = original_outcome)
    
    rocobj <- roc(outcome ~ predicted_num, 
                  data = prediction,
                  direction = '<',
                  levels = c(0, 1))
    
    fitted <- data.frame(
      mtry = optim_mtry,  
      splitrule = optim_splitrule,
      min.node.size = optim_min.node.size,
      auc = rocobj$auc,
      accuracy = mean(prediction$predicted==prediction$outcome),
      accuracy_num=mean(ifelse(prediction$predicted_num > 0.5, 1, 0) == prediction$outcome)
      )
    
    # BOOTSTRAP
    sampled_data <- sampler(data, 
                            outcome = outcome,
                            clust_var = clust_var,
                            sample_method = sample_method,
                            N = N, 
                            seed = seed)
    
    if(!is.null(clust_var)) clust_var <- "id"
    
    for (i in 1:N){
      
      ## sampled data in a matrix form
      sampled_outcome <- as.matrix(sampled_data[[1]][[i]]$outcome)
      
      if (!is.numeric(sampled_outcome)){
        sampled_outcome <- factor(sampled_outcome)
        mapping <- setNames(levels(sampled_outcome), 0:(length(levels(sampled_outcome))-1))
        sampled_outcome <- as.numeric(factor(sampled_outcome)) - 1 
      }
      
      sampled_outcome <- as.factor(sampled_outcome)
      
      sampled_predictors <- sampled_data[[1]][[i]] %>% 
        dplyr::select(-dplyr::all_of(c("outcome","id", "obs_id","Country", clust_var))) %>%
        as.matrix()
      
      ## re-optimize parameters
      # OPTIMALIZATION
      
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      
      rf_model_sampled <- caret::train(
        x = sampled_predictors, 
        y = sampled_outcome, 
        method = "ranger",  # Use kknn method for more flexibility
        tuneGrid = tune_grid,  # Manhattan or Euclidean distance
        metric = "Accuracy", # Specify the evaluation metric
        trControl = fitControl,
        importance = "impurity"
      )
      
      stopCluster(cl)
      
      optim_mtry_sampled = rf_model_sampled$bestTune$mtry
      optim_splitrule_sampled = rf_model_sampled$bestTune$splitrule
      optim_min.node.size_sampled = rf_model_sampled$bestTune$min.node.size
      
      # Combine predictors and outcome into a single data frame
      train_data_sampled <- data.frame(sampled_outcome, 
                                       sampled_predictors, 
                                       check.names = TRUE)
      # TRAINING
      sampled_fit <- ranger(
        formula = sampled_outcome ~ .,    # Explicitly specify the outcome variable 'y'
        data = train_data_sampled,        # Use the combined data frame
        mtry = optim_mtry_sampled,  # Use the best value of k from tuning
        splitrule = optim_splitrule_sampled,
        min.node.size = optim_min.node.size_sampled,
        probability = TRUE
      )
      
      # VALIDATION
      valid_outcome <- as.factor(as.matrix(sampled_data[[2]][[i]]$outcome))
      valid_predictors <- sampled_data[[2]][[i]] %>% 
        dplyr::select(-outcome, -id, - Country, -obs_id) %>% 
        as.matrix()
      
      valid_predicted_num = predict(sampled_fit, 
                                    data = as.data.frame(valid_predictors))$predictions[,2]
      
      valid_predicted_orig <- ifelse(valid_predicted_num > 0.5, 1, 0)
      
      
      prediction_onValidation <- data.frame(
        predicted = valid_predicted_orig,
        predicted_num = valid_predicted_num,
        outcome = valid_outcome,
        iteration = i)
      
      ## record performances
      kfold_rocobjs[[i]] <- roc(outcome ~ predicted_num, 
                                data = prediction_onValidation,
                                direction = '<',
                                levels = c(0, 1))
      
      auc_validation[i] <- kfold_rocobjs[[i]]$auc
      accuracy_validation[i] <- mean(prediction_onValidation$predicted == prediction_onValidation$outcome)
      predictions[[i]] <- data.frame(prediction_onValidation, id = sampled_data[[2]][[i]]$id)
    }
    
    ## connect predictions
    predictions <- bind_rows(predictions)
    
    ## aggregate obtained information
    valid_performances <- data.frame(
      auc_validation,
      accuracy_validation)
    
    model_summary <- fitted %>% 
      mutate(
        auc_optimism_corrected = mean(valid_performances$auc_validation),
        auc_optimism_corrected_CIL = quantile(valid_performances$auc_validation, probs = 0.025),
        auc_optimism_corrected_CIU = quantile(valid_performances$auc_validation, probs = 0.975),
        accuracy_optimism_corrected = mean(valid_performances$accuracy_validation),
        accuracy_optimism_corrected_CIL =  quantile(valid_performances$accuracy_validation, probs = 0.025),
        accuracy_optimism_corrected_CIU =  quantile(valid_performances$accuracy_validation, probs = 0.975)) %>% 
      dplyr::select(mtry,splitrule, min.node.size, auc,  
                    auc_optimism_corrected:auc_optimism_corrected_CIU,
                    accuracy,
                    accuracy_optimism_corrected:accuracy_optimism_corrected_CIU)
    
    ## define outputs
    rf_model <- list(model_summary = model_summary, 
                      valid_performances = valid_performances, 
                      predictions = prediction, 
                      rocobj=rocobj,
                      kfold_rocobjs=kfold_rocobjs,
                      trained_model=fit)
    
    # save results
    if (overfitting_check){
     if (!dir.exists(file.path("../intermediate_files/models_overfitting_check/",Q,file))){
      dir.create(file.path("../intermediate_files/models_overfitting_check/",Q,file))
     } 
      if (ra) save(rf_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model_ra.RData"))
      else save(rf_model,file=file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model.RData"))
    } else {
      if (!dir.exists(file.path("../intermediate_files/models/",Q,file))){
        dir.create(file.path("../intermediate_files/models/",Q,file))
      }
      if (ra) save(rf_model,file=file.path("../intermediate_files/models/",Q,file,"rf_model_ra.RData"))
      else save(rf_model,file=file.path("../intermediate_files/models/",Q,file,"rf_model.RData"))
    }
  }
  
  return(rf_model)
}

# other --------------
umap_plot_other <- function(asv_table,metadata, distance="robust.aitchison",neighbors=10,
                      min_dist=0.1,by=NULL){
  custom.config <- umap.defaults
  custom.config$random_state <- 123
  custom.config$n_neighbors <- neighbors
  custom.config$min_dist <- min_dist
  
  r.dist <- as.matrix(vegdist(t(asv_table %>% column_to_rownames("SeqID")), method=distance))
  umap_plot <- umap(r.dist,input="dist",config=custom.config)
  data_umap <- data.frame(umap_plot$layout)
  data_umap <- merge(data_umap %>% rownames_to_column("SampleID"),metadata,by="SampleID")
  if ("Country" %in% colnames(data_umap)){
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=Platform)) + geom_point() + theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") + 
      scale_color_manual(values=c("#f17d72", "#3bb4c4"))
    
  } else {
    p <- ggplot(data=data_umap, aes(x=X1,y=X2,color=Platform)) + geom_point() + theme_bw() +
      xlab("UMAP 1") + ylab("UMAP 2") + 
      scale_color_manual(values=c("#f17d72", "#3bb4c4"))
  }
  
  return(p)
}

pca_plots_OTHER <- function(mpse_object,by="Group"){
  plot_12 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    ellipse = TRUE, show.legend = FALSE ) +
    scale_fill_manual(
      values=c("#309f87","#425387"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#309f87","#425387"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC1 vs PC2")
  
  plot_13 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    .dim=c(1,3),
    ellipse = TRUE, show.legend = FALSE ) +
    scale_fill_manual(
      values=c("#309f87","#425387"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#309f87","#425387"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC1 vs PC3")
  
  plot_23 <- mpse_object %>% mp_plot_ord(
    .ord = pcoa, 
    .group = !!sym(by), 
    .color = !!sym(by), 
    .dim=c(2,3),
    ellipse = TRUE, show.legend = FALSE ) +
    scale_fill_manual(
      values=c("#309f87","#425387"), 
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_color_manual(
      values=c("#309f87","#425387"),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) +
    scale_size_continuous(
      range=c(0.5, 3),
      guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))) + 
    ggtitle("PC2 vs PC3")
  
  return(list(plot_12,plot_13,plot_23))
  
}

alpha_div_plot <- function(alpha_data){
  p_values <- c()
  for (index in c("Observe","Shannon","Simpson","Pielou")){
    res <- wilcox.test(alpha_data[,index][alpha_data$Platform=="ILLUMINA"],
                       alpha_data[,index][alpha_data$Platform=="AVITI"],
                       paired=TRUE) 
    p_values <- c(p_values,res$p.value)
  }
  
  sig <- p_values
  
  sig[p_values > 0.05] <-'N.S.'
  sig[p_values <= 0.05] <- 'P<0.05'
  sig[p_values <= 0.01] <- 'P<0.01'
  sig[p_values <= 0.001] <- 'P<0.001'
  
  alpha_data_melted <- melt(alpha_data)
  
  graphLabels <- data.frame(variable = c("Observe","Shannon","Simpson","Pielou"), Pval = sig)
  
  p <- ggplot() +
    geom_line(position=position_jitter(w=0.1, h=0,seed = 123),data=alpha_data_melted,aes(x = Platform,y = value,group = Pair),
              color = "grey",size = 0.6, alpha = 0.3) + 
    geom_boxplot(data=alpha_data_melted, aes(x=Platform, y=value),outliers = FALSE) + 
    geom_jitter(data=alpha_data_melted, aes(x=Platform, y=value, color=Platform),position = position_jitter(w=0.1, h=0,seed = 123)) + 
    scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) +
    scale_fill_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) + 
    # scale_fill_manual(values=c("#f9847c","#1ac6ca")) + 
    # scale_color_manual(values=c("#f9847c","#1ac6ca")) + 
    guides(fill="none",color = "none") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    facet_wrap(~factor(variable,levels=c("Observe","Shannon","Simpson","Pielou")), ncol = 4,scales = "free")  + 
    
    geom_text(data = graphLabels, aes(Inf, Inf, label = Pval),hjust = 1.2, vjust = 1.2,size=5, show.legend = FALSE) 
  
  return(p)
}


beta_div_plot <- function(ps,metadata,rarefy=TRUE,normalize=FALSE,filter=FALSE,measure="aitchison"){
  if(rarefy){
    ps <- rarefy_even_depth(ps,sample.size = 10000)
  }
  if (filter){
    filt_data <- abundance_filtering(abundance_threshold = 0.01,
                                     as.data.frame(ps@otu_table) %>% rownames_to_column("SeqID"),
                                     as.data.frame(ps@tax_table) %>% rownames_to_column("SeqID"))
    
    ps <- construct_phyloseq(filt_data[[1]],filt_data[[2]],metadata)
  }
  if(normalize){
    ps <- transform_sample_counts(ps, function(x) x/sum(x))
  }
  if (measure=="aitchison") {
    ps <- microbiome::transform(ps,"clr")
    measure_final <- "euclidean"
  } else if (measure=="robust.aitchison"){
    ps <- microbiome::transform(ps,"rclr")
    measure_final <- "euclidean"
  } else measure_final <- measure
  
  ord <- ordinate(ps, method = "PCoA", distance = measure_final)
  imp_vec <- ord$values$Relative_eig
  pca_vec <- ord$vectors
  lab1 = paste("PCo1 ", "(",round(imp_vec[1]*100,2),"%", ")", sep="")
  lab2 = paste("PCo2 ", "(",round(imp_vec[2]*100,2),"%", ")", sep="")
  lab3 = paste("PCo3 ", "(",round(imp_vec[3]*100,2),"%", ")", sep="")
  
  data_for_pca <- as.data.frame(t(asv_tab[,-which(colnames(asv_tab)=="SeqID")]))
  data_for_pca <- cbind(data_for_pca,metadata)
  
  p_12 <- ggplot(data_for_pca, aes(x=pca_vec[,1],y=pca_vec[,2],color=!!sym("Platform"))) + 
    geom_line(aes(group=Pair), colour="grey") +
    geom_point(show.legend = TRUE,size=3) + stat_ellipse(show.legend = FALSE) + 
    xlab(lab1)+
    ylab(lab2)+
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(legend.position = "right") + 
    scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) +
    scale_fill_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) 
  #scale_fill_manual(values=c("#f9847c","#1ac6ca")) + 
  #scale_color_manual(values=c("#f9847c","#1ac6ca")) 
  
  p_13 <- ggplot(data_for_pca, aes(x=pca_vec[,1],y=pca_vec[,3],color=!!sym("Platform"))) + 
    geom_line(aes(group=Pair), colour="grey") +
    geom_point(show.legend = TRUE, size=3) + stat_ellipse(show.legend = FALSE) + 
    xlab(lab1)+
    ylab(lab3)+
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) +  
    theme(legend.position = "right") + 
    #  scale_fill_manual(values=c("#f9847c","#1ac6ca")) + 
    #   scale_color_manual(values=c("#f9847c","#1ac6ca")) 
    scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) +
    scale_fill_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) 
  
  p_23 <- ggplot(data_for_pca, aes(x=pca_vec[,2],y=pca_vec[,3],color=!!sym("Platform"))) + 
    geom_line(aes(group=Pair), colour="grey") +
    geom_point(show.legend = TRUE,size=3) + stat_ellipse(show.legend = FALSE) + 
    xlab(lab2)+
    ylab(lab3)+
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(legend.position = "right") + 
    # scale_fill_manual(values=c("#f9847c","#1ac6ca")) + 
    #   scale_color_manual(values=c("#f9847c","#1ac6ca")) 
    scale_color_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) +
    scale_fill_manual(values=c("#309f87","#425387","#f9c675","#d55c4a")) 
  
  p <- ggarrange(p_12,p_13,p_23,ncol=3)
  return(p_12)
}



mock_zymo_genus_merging <- function(asv_table,taxa_table,zymo_asv_table, zymo_taxa){
  # merging mock community samples with reference abundance
  taxa_reads_table <- create_asv_taxa_table(asv_table,taxa_table)
  genus_data <- aggregate_taxa(asv_table,taxa_table, taxonomic_level = "Genus")
  genus_asv_table <- genus_data[[1]]
  genus_taxa_table <- genus_data[[2]]
  
  genus_asv_table_norm <- as.data.frame(apply(genus_asv_table[,-1],2,function(x) x/sum(x)))
  genus_asv_table_norm$SeqID <- genus_asv_table$SeqID
  
  colnames(zymo_asv_table)[1] <- "SeqID"
  colnames(zymo_taxa)[1] <- "SeqID"
  
  genus_data_zymo <- aggregate_taxa(zymo_asv_table,zymo_taxa, taxonomic_level = "Genus")
  genus_zymo_asv_table <- genus_data_zymo[[1]][-9,]
  genus_zymo_taxa_table <- genus_data_zymo[[2]][-9,]
  
  merged_data <- merge(genus_zymo_asv_table,genus_asv_table_norm, by="SeqID", all=TRUE)
  merged_data[is.na(merged_data)] <- 0
  colnames(merged_data)[2] <- "ZYMO REFERENCE"
  return(merged_data)
}


precision_recall <- function(mock_zymo_genus){
  # calculating precision
  precs <- c()
  recs <- c()
  zymo_taxa_df <- mock_zymo_genus$`ZYMO REFERENCE`
  zymo_taxons <- zymo_taxa_df>0
  zymo_taxons <- mock_zymo_genus$SeqID[zymo_taxons]
  for (sample in colnames(mock_zymo_genus)[-c(1,2)]){
    sample_taxa_df <- mock_zymo_genus[,sample]
    sample_taxons <- sample_taxa_df>0
    sample_taxons <- mock_zymo_genus$SeqID[sample_taxons]
    
    tp <- sum(sample_taxons %in% zymo_taxons)
    fp <- sum(!(sample_taxons %in% zymo_taxons))
    fn <- sum(!(zymo_taxons %in% sample_taxons))
    rec = tp/(tp+fn)
    recs <- c(recs,rec)
    
    prec = tp/(tp + fp)
    precs <- c(precs,prec)
  }
  names(precs) <- colnames(mock_zymo_genus)[-c(1,2)]
  names(recs) <- colnames(mock_zymo_genus)[-c(1,2)]
  df <- data.frame(Precision=precs,
                   Recall=recs)
  return(df)
}


composition <- function(ps, rarefy=FALSE, legend=FALSE, 
                        legend_position="right",groupby=NULL){
  if (rarefy) {
    ps <- rarefy_even_depth(ps,sample.size = 10000)
  }
  ps <- suppressWarnings(merge_samples(ps, groupby))
  ps <- transform_sample_counts(ps, function(x) x/sum(x))
  
  p <- plot_bar(ps, fill = "Phylum") +
    geom_bar(aes(color=Phylum, 
                 fill=Phylum), stat="identity", position="stack") + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    scale_color_manual(breaks=c("Actidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibionota","Campilobacterota","Chloroflexi","Cyanobacteria","Deinococcota","Dependentiae","Desulfobacterota","Elusimicrobiota","Euryarcheota","Firmicutes","Fusobacteriota","Myxzcoccota","NB1-j","Patescibacteria","Planctomycetota","Proteobacteria","RCP2-54","Spirochaetota","Synergistota","Thermoplasmatota","unassigned","Verrucomicrobiota","WPS-2"),
                       values=c("#1f1f1f","#546494","#388842","#8350b6","#d62e5f","#381216","#6e94c8","#aeae43","#9d729f","#ba3237","#723f23","#4f9eba","#b3dd9c","#e18dac","#ae4a37","#dbad80","#8ebcca","#e7d836","#e7928f","#dc6438","#d7c2a3","#e9e07f","#d1a9b0","#e8b14f","#309f87","#425387","#f9c675","#d55c4a"))+
    
    
    scale_fill_manual(breaks=c("Actidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibionota","Campilobacterota","Chloroflexi","Cyanobacteria","Deinococcota","Dependentiae","Desulfobacterota","Elusimicrobiota","Euryarcheota","Firmicutes","Fusobacteriota","Myxzcoccota","NB1-j","Patescibacteria","Planctomycetota","Proteobacteria","RCP2-54","Spirochaetota","Synergistota","Thermoplasmatota","unassigned","Verrucomicrobiota","WPS-2"),
                      values=c("#1f1f1f","#546494","#388842","#8350b6","#d62e5f","#381216","#6e94c8","#aeae43","#9d729f","#ba3237","#723f23","#4f9eba","#b3dd9c","#e18dac","#ae4a37","#dbad80","#8ebcca","#e7d836","#e7928f","#dc6438","#d7c2a3","#e9e07f","#d1a9b0","#e8b14f","#309f87","#425387","#f9c675","#d55c4a"))
  
  if (!legend) p <- p + theme(legend.position = "none")
  else p <- p + theme(legend.position = legend_position)
  return (p)
}

mock_zymo_genus_plot <- function(mock_zymo_genus, metadata, setting){
  
  metadata <- metadata %>% filter(Sample %in% colnames(mock_zymo_genus)[-1]) %>% 
    column_to_rownames("Sample")
  
  mock_otf <- rownames(metadata[metadata$primers=="OTF",])
  mock_stf <- rownames(metadata[metadata$primers=="STF",])
  
  mock_zymo_genus <- mock_zymo_genus[,c("SeqID", "ZYMO REFERENCE", mock_stf,mock_otf)]
  
  
  tree <- hclust(dist(t(mock_zymo_genus[,-1])))
  tree <- dendro_data(tree)
  
  metadata <- metadata[colnames(mock_zymo_genus)[-1],c("primers","sample_type","run_number",setting)]
  metadata[1,] <- c("ZR","ZR","ZR","ZR")
  rownames(metadata)[1] <- "ZYMO REFERENCE"
  
  mock_zymo_genus_samples <- mock_zymo_genus
  metadata_samples <- metadata
  mock_zymo_genus_melted_samples <- melt(mock_zymo_genus_samples)
  colnames(mock_zymo_genus_melted_samples) <- c("SeqID", "Sample", "Relative abundance")
  
  mock_zymo_genus <- mock_zymo_genus[,c("SeqID",tree$labels$label)]
  metadata <- metadata[colnames(mock_zymo_genus)[-1],c("primers","sample_type","run_number",setting)]
  
  mock_zymo_genus_melted <- melt(mock_zymo_genus)
  colnames(mock_zymo_genus_melted) <- c("SeqID", "Sample", "Relative abundance")
  
  
  if (!(FALSE %in% (rownames(metadata) == colnames(mock_zymo_genus)[-1]))){
    labels <- paste(metadata$primers,metadata$run_number,as.character(round(as.numeric(metadata[[setting]]),2)), sep="\n")
  labels_samples <- paste(metadata_samples$primers,metadata_samples$run_number,as.character(round(as.numeric(metadata_samples[[setting]]),2)), sep="\n")
    } 
  
  p <- ggplot() + 
    geom_bar(data=mock_zymo_genus_melted, aes(y=`Relative abundance`, x=Sample, fill=SeqID),position="fill", stat="identity") + 
    geom_text(size=2,aes(label = labels,x = 1:length(labels),y=rep(1.07,length(labels)))) + 
    scale_y_continuous(limits = c(-0.2,1.1),breaks=seq(0,1,by=0.2)) + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(legend.position = "right",axis.text.x = element_blank(),
          legend.text=element_text(size=10),
          legend.title = element_text(size=15),  
          axis.title.y = element_text(size=15), 
          axis.title.x = element_text(size=15),
          axis.text.y=element_text(size=10),
          axis.ticks = element_blank())  + 
    scale_fill_manual(values=setNames(c("#546494","#388842","#8350b6","#d62e5f","#381216","#6e94c8","#aeae43","#9d729f"),c("Pseudomonas","Escherichia-Shigella","Salmonela","Lactobacillus","Enterococcus","Staphylococcus","Listeria","Bacillus"))) + 
    guides(fill=guide_legend(title="Genus")) + 
    geom_segment(
      data = tree$segments,
      aes(x = x, y = -y, xend = xend, yend = -yend)
    )
  
  p_samples <- ggplot() + 
    geom_bar(data=mock_zymo_genus_melted_samples, aes(y=`Relative abundance`, x=Sample, fill=SeqID),position="fill", stat="identity") + 
    geom_text(size=2,aes(label = labels_samples,x = 1:length(labels_samples),y=rep(1.07,length(labels_samples)))) + 
    scale_y_continuous(limits = c(0,1.1),breaks=seq(0,1,by=0.2)) + 
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0)) + 
    theme(legend.position = "right",axis.text.x = element_blank(),
          legend.text=element_text(size=10),
          legend.title = element_text(size=15),  
          axis.title.y = element_text(size=15), 
          axis.title.x = element_text(size=15),
          axis.text.y=element_text(size=10),
          axis.ticks = element_blank())  + 
    scale_fill_manual(values=setNames(c("#546494","#388842","#8350b6","#d62e5f","#381216","#6e94c8","#aeae43","#9d729f"),c("Pseudomonas","Escherichia-Shigella","Salmonela","Lactobacillus","Enterococcus","Staphylococcus","Listeria","Bacillus"))) + 
    guides(fill=guide_legend(title="Genus")) 
    
  
  
  return(list(p,p_samples))
  
}
