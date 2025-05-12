suppressMessages(suppressWarnings({
  library(data.table)
  library(ccrepe)
  library(igraph)
  library(readr)
  library(cowplot)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(phyloseq)
  library(MicrobiotaProcess)
  library(ggpubr)
  library(ggrepel)
  library(ggplotify)
  #library(radEmu)
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
  #library(umap)
  library(Maaslin2)
  library(ggvenn)
  library(ranger)
  library(doParallel)
  library(gbm)
  library(tidyr)
  library(kableExtra)
  library(tidyverse)
  library(picante)
  library(tidyr)
  library(mice)
}))

## Basic functions for data processing ----

read_counts <- function(asv_table, text=FALSE, line=5000){
  # Creates a plot of library size - read count per sample
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
  # Checks the data from redundant taxa
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # outputs:
  # list(asv_table, taxa_table) - non-redundant tables
  
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
  
  # Creates the dataset by merging two cohorts by specific segment
  # inputs:
  # asv_tab_1/2 - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_tab_1/2 - Taxonomy with 'SeqID' and taxonomy
  # metadata_1/2 - metadata with 'SampleID' identifier
  # segment - terminal_ileum/colon
  # Q - analysis question (Q1/Q2/Q3)
  # outputs:
  # list(asv_table, taxa_table, metadata) - newly created dataset
  
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
    if (Q=="Q5"){
     merged_metadata <- rbind(metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country, Calprotectin),
                             metadata_2 %>% dplyr::select(SampleID,subjectid,Group,segment,Country, Calprotectin) %>% 
                               dplyr::rename(Patient=subjectid,Matrix=segment) %>%
                               mutate(Patient=paste0("NO_",Patient))) 
    } else{
      merged_metadata <- rbind(metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country),
                               metadata_2 %>% dplyr::select(SampleID,subjectid,Group,segment,Country) %>% dplyr::rename(Patient=subjectid,Matrix=segment) %>%
                                 mutate(Patient=paste0("NO_",Patient))) 
    }
    
    row.names(merged_metadata) <- NULL
    
  } else {
    # Merging
    merged_asv_tab <- asv_tab_1
    merged_taxa_tab <- taxa_tab_1
    
    # metadata merge
    if ("Patient" %in% colnames(metadata_1)){
      if (Q=="Q5")  merged_metadata <- metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country, Calprotectin)
      else merged_metadata <- metadata_1 %>% dplyr::select(SampleID,Patient,Group,Matrix,Country)
    } else {
      if (Q=="Q5") {merged_metadata <- metadata_1 %>% 
        dplyr::select(SampleID,subjectid,Group,Matrix,Country, Calprotectin) %>% 
        dplyr::rename(Patient=subjectid) %>%
        mutate(Patient=paste0("NO_",Patient))}
      else {
        merged_metadata <- metadata_1 %>% 
          dplyr::select(SampleID,subjectid,Group,Matrix,Country) %>% 
          dplyr::rename(Patient=subjectid) %>%
          mutate(Patient=paste0("NO_",Patient))
      }
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
    merged_metadata <- merged_metadata[merged_metadata$Group %in% c("rPSC","non-rPSC", "pre_ltx"),]
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
  # Creates the dataset by merging two cohorts by specific segment
  # inputs:
  # asv_tab_1/2 - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_tab_1/2 - Taxonomy with 'SeqID' and taxonomy
  # metadata_1/2 - metadata with 'SampleID' identifier
  # segment - terminal_ileum/colon
  # Q - analysis question (Q1/Q2/Q3)
  # outputs:
  # list(asv_table, taxa_table, metadata) - newly created dataset
  
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
  # Creates the dataset by merging two cohorts by specific segment
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # outputs:
  # taxa_asv_table - data frame with counts and taxonomy
  
  # creating asv+taxa
  asvs <- asv_table$SeqID
  taxa_ranks <- colnames(taxa_table)
  where_level <- which(tolower(taxa_ranks)=="species")
  if (length(where_level)==0) where_level <- which(tolower(taxa_ranks)=="genus")
  if (length(where_level)==0) where_level <- which(tolower(taxa_ranks)=="phylum")
  if (length(where_level)==0) where_level <- which(tolower(taxa_ranks)=="domain")
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  
  if (where_level!=2){
    seq_ids <- apply(taxa_asv_table[,2:where_level],1, function(x){
    a <- paste0(substring(tolower(colnames(taxa_asv_table[,2:where_level])),1,1),"__",x, collapse = ";")
    return(a)
  })
  } else {
    seq_ids <- paste0("gu__",taxa_asv_table[,2])
  }
  
  
  taxa_asv_table$Taxonomy <- seq_ids
  taxa_asv_table %<>% column_to_rownames("SeqID")
  taxa_asv_table <- taxa_asv_table[asvs,]
  taxa_asv_table <- taxa_asv_table[,-which(colnames(taxa_asv_table) %in% taxa_ranks[2:8])]
  taxa_asv_table <- taxa_asv_table[,c(ncol(taxa_asv_table),1:(ncol(taxa_asv_table)-1))]
  colnames(taxa_asv_table) <- c("SeqID", colnames(taxa_asv_table)[-1])
  return(taxa_asv_table)
}


binomial_prep <- function(asv_table,taxa_table,metadata,group, patient=FALSE,
                          usage="linDA",filtering=TRUE){
  # Prepares the dataset for univariate analysis or machine learning functions
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # metadata - data frame with 'SampleID' identifier
  # group - what groups will be tested, e.g. (c(healthy,pre_LTx))
  # patient - boolean, should info about patient be incorporated?
  # usage - linDA/ml_clr/ml_ra 
  # outputs:
  # list(uni_data,uni_tax,uni_metadata) or 
  # uni_data
  
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
  
  if (filtering){
    filt_data <- filtering_steps(uni_data %>% rownames_to_column("SeqID"),uni_tax,
                               uni_metadata %>% rownames_to_column("SampleID"),
                               seq_depth_threshold=10000)
  uni_data <- filt_data[[1]]  %>% column_to_rownames("SeqID")
  uni_tax <- filt_data[[2]]
  uni_metadata <- filt_data[[3]] %>% column_to_rownames("SampleID")
  }  else{
    data_filt <- seq_depth_filtering(uni_data %>% rownames_to_column("SeqID"),uni_tax,
                                     uni_metadata %>% rownames_to_column("SampleID"),
                                     seq_depth_threshold=10000)
    uni_data <- data_filt[[1]]  %>% column_to_rownames("SeqID")
    uni_tax <- data_filt[[2]]
    uni_metadata <- data_filt[[3]] %>% column_to_rownames("SampleID")
  }
  
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
  # This functions is similar to binomial_prep(), 
  # but works for the PSC effect data preparation
  
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


## Filtering functions ----

filtering_steps <- function(asv_tab,taxa_tab,metadata,
                            seq_depth_threshold=10000){
  
  # Filters the dataset by predefined rules 
  # (sequencing_depth, nearzerovar()) and checks the low depth of samples
  # afterwards
  # inputs:
  # asv_tab - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_tab - Taxonomy with 'SeqID' and taxonomy
  # metadata -  metadata with 'SampleID' identifier
  # seq_depth_threshold - threshold for sequencing depth filtering, default 10,000
  # outputs:
  # list(filt_asv_tab,filt_taxa_tab, filt_metadata)
  
  # filtering
  ## sequencing depth
  data_filt <- seq_depth_filtering(asv_tab,taxa_tab,metadata,
                                   seq_depth_threshold = seq_depth_threshold)
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
  # Performs the sequencing depth filtering
  # inputs:
  # asv_tab - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_tab - Taxonomy with 'SeqID' and taxonomy
  # metadata -  metadata with 'SampleID' identifier
  # seq_depth_threshold - threshold for sequencing depth filtering, default 10,000
  # outputs:
  # list(filt_asv_tab,filt_taxa_tab, filt_metadata)
  
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

nearzerovar_filtering <- function(asv_tab,taxa_tab,metadata){
  # Performs the nearzerovar() filtering
  # inputs:
  # asv_tab - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_tab - Taxonomy with 'SeqID' and taxonomy
  # metadata -  metadata with 'SampleID' identifier
  # outputs:
  # list(filt_asv_tab,filt_taxa_tab)
  
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
  # Checks the dataset for low depth samples and removes them if < 100
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # metadata -  metadata with 'SampleID' identifier
  # outputs:
  # list(filt_asv_table, filt_taxa_tab, filt_metadata)
  
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
  # Creates the summary of filtering steps
  # inputs:
  # asv_tab - original data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # filt_asv_tab - filtered data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # filt_metadata - filtered metadata with 'SampleID' identifier
  # seq_step - number of taxa in filtered dataset
  # prev_step - 0
  # nearzero_step - number of taxa in filtered dataset (after nearzerovar())
  # outputs:
  # df - dataframe with summarized info about number of ASVs, samples, 
  # number of filtered ASVs, samples and so on
  
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

pairwise.lm <- function(formula,factors,data, p.adjust.m ='BH')
{
  # Runs linear model (LM) for each pairwise comparison of groups in input and
  # performs the BH correction (or any other set by user). 
  # In this analysis, it is used in alpha diversity testing. That is why 
  # the alpha diversity indices are required as input
  # inputs:
  # formula - formula for LM
  # factors - vector or groups to be tested (column from data dataframe)
  # data - data frame with computed alpha diversity indices
  # p.adjust.m - method for p adjustment, default BH
  # outputs:
  # (list(df,emeans_models, means)): dataframe with results, emeans results when 
  # posthoc analysis was needed
  
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
  if (all(c("healthy","non-rPSC","rPSC","pre_ltx") %in% factors)){
    co <- data.frame("1"=c("healthy","healthy","healthy","non-rPSC","pre_ltx","pre_ltx"),
                     "2"=c("non-rPSC","rPSC","pre_ltx","rPSC","non-rPSC","rPSC")) %>% t()
  }
  if ("ibd" %in% factors) {
    co <- data.frame("1"=c("no_ibd"),
                     "2"=c("ibd")) %>% t()
  } 
  if ("high" %in% factors) {
    co <- data.frame("1"=c("low"),
                     "2"=c("high")) %>% t()
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
  if (!is.data.frame(models)) models %<>% as.data.frame() %>% `row.names<-`(names[-grep("^[(].+[)]$", names)])
  if (ncol(co)>1){
    p.adjusted <- p.adjust(models[,"Pr(>|t|)"],method=p.adjust.m)
    models$p.adj <- p.adjusted
    sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  df <- data.frame(models,
                   sig=sig)
  } else df <- models
  
  
  return(list(df,emeans_models, means))
  
}

pairwise.lmer <- function(formula,factors,data, p.adjust.m ='BH')
{
  # This function is similar to pairwise.lm(), but performs linear mixed effect models
  # instead of simple linear model to include the random effect of patient in this analysis
  
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
  if ("low" %in% factors) {
    co <- data.frame("1"=c("low"),
                     "2"=c("high")) %>% t()
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
      if (model[4,"Pr(>|t|)"]<0.1){
      # CZ
      group1_group2_cz <- as.data.frame(coef(summary(lmer(as.formula(gsub(" \\* Country","",formula)), data = subset(x_sub,Country=="CZ")))))
      rownames(group1_group2_cz)[1] <- paste0(co[1,elem],"_CZ")
      # NO
      group1_group2_no <- as.data.frame(coef(summary(lmer(as.formula(gsub(" \\* Country","",formula)), data = subset(x_sub,Country=="NO")))))
      rownames(group1_group2_no)[1] <- paste0(co[1,elem],"_NO")
      # GROUP 1 
      group1_country <- as.data.frame(coef(summary(lmer(as.formula(gsub(" Group \\*","",formula)), data = subset(x_sub,Group==co[1,elem])))))
      rownames(group1_country)[1] <- paste0(co[1,elem],"_CZ")
      # GROUP 2 
      group2_country <- as.data.frame(coef(summary(lmer(as.formula(gsub(" Group \\*","",formula)), data = subset(x_sub,Group==co[2,elem])))))
      rownames(group2_country)[1] <- paste0(co[2,elem],"_CZ")
      
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

glm_renaming <- function(glm_data, group){
  # Renames the results of linDA (linda() function $output) based on the chosen group[1] 
  # inputs:
  # linda_data - result of linda()$output
  # group - strings corresponding to the name of the group[1]
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
  # Renames the results of linDA (linda() function $output) based on the chosen group[1] 
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

pairwise.adonis <- function(x,factors, covariate=NULL, interaction=FALSE,patients=NULL, 
                            sim.function = 'vegdist', sim.method = 'robust.aitchison',
                            p.adjust.m ='BH',perm=999)
{
  # Runs adonis() function (PERMANOVA) for each pairwise comparison of groups in input and
  # performs the BH correction (or any other set by user). 
  # Main effects of Group and Country are tested using the by = 'margin' setting. 
  # Interaction effect is calculated with by = 'terms' setting.
  # inputs:
  # x - asv table with SeqID as colnames and samples as rownames
  # factors - vector or groups to be tested (column from data dataframe)
  # covariate - vector of covariates (for example Country) that needs to be included as main effect
  # interaction - boolean, if interaction effect should be calculated, default FALSE
  # patients - vector of patients, which will be used for custom permutations, default NULL
  # sim.function - function to be used for distance calculation, default vegdist
  # sim.method - distance metric, default robust.aitchison
  # perm - number of permutations, if Patients are provided, custom permutations will be generated,
  # p.adjust.m - method for p adjustment, default BH
  # outputs:
  # list(pairw.res_factor,pairw.res_cov,pairw.res_fac.cov): 
  # results for main effects of factor, covariate and interaction effect
  
  set.seed(123)
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
  if ("high" %in% factors) {
    co <- data.frame("1"=c("low"),
                     "2"=c("high")) %>% t()
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
  # Performs post-analysis of pairwise PERMANOVA.
  # In terms of PSC_study, it runs PERMANOVA individually for groups comparison 
  # for each country individually
  # inputs:
  # x - asv table with SeqID as colnames and samples as rownames
  # factors - vector or groups to be tested (column from data dataframe)
  # covariate - vector of covariates (for example Country) that needs to be included as main effect
  # group1 - name of the first group
  # group2 - name of the second group
  # patients - vector of patients, which will be used for custom permutations, default NULL
  # sim.method - distance metric, default robust.aitchison
  # outputs:
  # list(group1_group2_cz,group1_group2_no,group1_country,group2_country)
  
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

custom_permutations <- function(x,factors,patients, perm=999){
  # Generates custom permutations, used in adonis()
  
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

## Visualization functions ----

taxon_boxplot <- function(asv_tab,taxa_tab,metadata,taxon){
  asv_tab <- as.data.frame(rrarefy(
    asv_tab %>% t(), 
    sample=10000)) %>% t() %>% as.data.frame() %>% rownames_to_column("SeqID")
  
  boxplot_df <- asv_tab %>% dplyr::filter(SeqID==taxon) %>% 
    column_to_rownames("SeqID") %>% t() %>% as.data.frame() %>%
    rownames_to_column("SampleID")
  
  boxplot_df <- merge(boxplot_df,metadata,by="SampleID",all.x=TRUE)
  boxplot_df %<>% dplyr::mutate(taxon=log(!!sym(taxon)))
  p <- ggplot(data=boxplot_df) + 
    geom_boxplot(aes(x=Group,y=(taxon)))
  
}

horizontal_barplot <- function(wb,taxa){
  # Generates a horizontal barplot showing the prevalence of taxa in individual groups,
  # This plot is visualized next to the dot-heatmap in DAA visualization
  # inputs:
  # wb - workbook that contains the univariate statistics about each taxon
  # taxa - taxa for which barplots should be generated, usually the one in dotheatmap
  # outputs:
  # p - horizontal barplot, ggplot object
  
  sheets_names <- sheets(wb)
  groups_names <- unique(unlist(strsplit(sheets_names,split=" vs ")))
  prevalences_df <- NULL
  
  for (name in sheets_names){
    group1_name <- unlist(strsplit(name,split=" vs "))[1]
    group2_name <- unlist(strsplit(name,split=" vs "))[2]
    df <- readWorkbook(wb, sheet = name)
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
  } else if (any(grepl("rPSC",colnames(prevalences_df))) & any(grepl("pre_ltx",colnames(prevalences_df)))){
    prevalences_df <- prevalences_df %>%
      group_by(SeqID) %>%
      summarise(
        PREVALENCE_.healthy = coalesce(PREVALENCE_.healthy[1], PREVALENCE_.healthy[2]),
        PREVALENCE_.pre_ltx = coalesce(PREVALENCE_.pre_ltx[1], PREVALENCE_.pre_ltx[2]),
        `PREVALENCE_.non-rPSC` = coalesce(`PREVALENCE_.non-rPSC`[1], `PREVALENCE_.non-rPSC`[2]),
        PREVALENCE_.rPSC = coalesce(PREVALENCE_.rPSC[1], PREVALENCE_.rPSC[2])
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
  else if (("rPSC" %in% prevalences_df_melt$variable) & ("pre_LTx" %in% prevalences_df_melt$variable)) {
    colors <- c("#309f87","#f9c675","#F08080","#A00000")
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("HC","pre_LTx", "non-rPSC","rPSC"))
  } else if ("rPSC" %in% prevalences_df_melt$variable) {
    colors <- c("#309f87","#F08080","#A00000")
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("HC","non-rPSC","rPSC"))
  } else if ("ibd" %in% alpha_data$Group){
    colors <- c("#A06A2C", "#B2182B")  
    
    prevalences_df_melt$variable <- factor(prevalences_df_melt$variable,levels = c("no_ibd","ibd"))
  }
  else (cat("chyba","\n"))
  
  p <- ggplot(prevalences_df_melt, aes(x = SeqID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    coord_flip() + 
    scale_fill_manual(values=colors)  +
    theme_minimal() + 
    theme(axis.title.x =element_text(size=8),
          axis.text.x = element_text(size=4),
          axis.text.y = element_blank(),
          axis.line.y.left = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.ticks.y.left = element_blank())+
    scale_y_continuous(breaks=c(0,0.5,1)) + theme(axis.text.x=element_text(size=8)) + 
    ylab("Prevalence") + 
    xlab("") + 
    theme(legend.position = "none")
  
  return(p)
}

alpha_diversity_custom_2 <- function(alpha_data, size=1.5,width=0.2){
  alpha_data %<>% drop_na()
  #colnames(alpha_data)[which(colnames(alpha_data)=="Observe")] <- "Richness"
  alpha_data <- alpha_data %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  
  if ("Observe" %in% colnames(alpha_data)) alpha_data %<>% dplyr::mutate(Richness=Observe)
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
  }
    else if ("high" %in% alpha_data$Group){
      colors <- c("#A06A2C", "#B2182B")  
      alpha_data$Group <- factor(alpha_data$Group,levels = c("low","high"))
  } else (print("chyba","\n"))
  richness_limit <- max(alpha_data$Richness) + 0.2*max(alpha_data$Richness)
  p_richness <- ggplot() + 
    geom_boxplot(data=alpha_data, aes(x=Group, y=Richness),outliers = FALSE) + 
    geom_jitter(width = width,height = 0,data=alpha_data,aes(x=Group, y=Richness, color=Group,
                                                             shape=Country),size=size) +
    scale_y_continuous(breaks = seq(0, richness_limit, by = 50)) + 
    ylim(0,richness_limit) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    xlab("") + 
    ylab("Richness") + 
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
  
  shannon_limit <- max(alpha_data$Shannon) + 0.2*max(alpha_data$Shannon)
  p_shannon <- ggplot() + 
    geom_boxplot(data=alpha_data, aes(x=Group, y=Shannon),outliers = FALSE) + 
    geom_jitter(width = width,height = 0,data=alpha_data,aes(x=Group, y=Shannon, color=Group,
                                                             shape=Country),size=size) +
    scale_y_continuous(breaks = seq(0, shannon_limit, by = 1)) + 
    ylim(0,shannon_limit) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    xlab("") + 
    ylab("Shannon") + 
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
  
  
  p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="none")
  return(p)
}
  
alpha_diversity_countries <-function(alpha_data,show_legend=FALSE){
  # Generates a scatterplot + boxplot showing the alpha diversity indices 
  # (Richness and Shannon) for each country individually
  # inputs:
  # alpha_data - metadata-like dataframe with 'SampleID' identifier
  #              and 'Observe' and 'Shannon' columns
  # show_legend - boolean, should legend be visualized?, default FALSE 
  #               (legends are added manually in publication-ready plots) 
  # outputs:
  # p - horizontal barplot, ggplot object
  
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
  } else if ("low" %in% alpha_data$Group){
    #colors <- c("#1B7837","#B2182B")  
    #colors <- c("#D04E36", "#3F7D3C")  
    colors <- c("#A06A2C", "#B2182B")  
    
    alpha_data$Group <- factor(alpha_data$Group,levels = c("low","high"))
  }
  
  else (print("chyba","\n"))
  
  alpha_data$Country <- factor(alpha_data$Country,levels = c('CZ','NO'))
  alpha_richness_data <- melt(alpha_data) %>% dplyr::filter(variable=="Observe")
  
  richness_limit <- max(alpha_data$Observe) + 0.2*max(alpha_data$Observe)
  
  p_richness <- ggplot(alpha_richness_data, aes(x=Group, y=value, fill=Country)) + 
    geom_boxplot(aes(fill=Country), outlier.shape = NA, position=position_dodge(width=0.8)) + 
    geom_jitter(aes(x=Group, y=value, color=Group, shape=Country), 
                position=position_jitterdodge(jitter.width=0.7, dodge.width=0.8), 
                size=0.5) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    xlab("") + 
    ylab("Richness") + 
    scale_y_continuous(breaks = seq(0, richness_limit, by = 50)) + 
    ylim(0,richness_limit)  + 
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
  
  alpha_shannon_data <- melt(alpha_data) %>% dplyr::filter(variable=="Shannon")
  shannon_limit <- max(alpha_data$Shannon) + 0.2*max(alpha_data$Shannon)
  
  p_shannon <- ggplot(alpha_shannon_data, aes(x=Group, y=value, fill=Country)) + 
    geom_boxplot(aes(fill=Country), outlier.shape = NA, position=position_dodge(width=0.8)) + 
    geom_jitter(aes(x=Group, y=value, color=Group, shape=Country), 
                position=position_jitterdodge(jitter.width=0.7, dodge.width=0.8), 
                size=0.5) + 
    scale_color_manual(values=colors) +
    scale_fill_manual(values=c("#ffffff","#ffffff")) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    xlab("") + 
    ylab("Shannon") + 
    scale_y_continuous(breaks = seq(0, shannon_limit, by = 50)) + 
    ylim(0,shannon_limit) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
  
  if (show_legend){
    p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="right")
  } else {
    p <- ggarrange(p_richness,p_shannon,common.legend = TRUE,legend="none")
  }
  
  return(p)
}

pca_plot_custom <- function(asv_table,taxa_table,metadata, 
                            measure="robust.aitchison",
                            show_boxplots = TRUE,
                            variable = "Group", size=2,
                            show_legend=TRUE,
                            clinical=FALSE,
                            clinical_metadata=NULL){
  # Generates a PCA plot showing the beta diversity of the dataset 
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # metadata - data frame with 'SampleID' identifier
  # measure - beta diversity distance metric, default robust.aitchison
  # show_boxplots - boolean, should boxplots be visualized along the PC, default TRUE 
  # show_legend - boolean, should legend be visualized?, default TRUE 
  # variable - variable that should be used for color-coding, default 'Group'
  # size - size of the points, default 2
  # clinical - boolean, if clinical variables should be visualized as BIPLOT, default FALSE
  # clinical_metadata, if clinical=TRUE, clinical_metadata should be provided with 'SampleID' identifier
  # outputs:
  # p - PCA plot, ggplot object
  
  if (variable=="Group"){
    metadata <- metadata %>%
    mutate(Group = case_when(
      Group == "healthy" ~ "HC",
      Group == "post_ltx" ~ "post_LTx",
      Group == "pre_ltx" ~ "pre_LTx",
      TRUE ~ Group  # Keep other values as is
    ))
  }
  
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
      dplyr::select(-c(Group,Country))
    
    #if (length(unique(clinical_metadata$PatientID))==nrow(clinical_metadata)){
    #  clinical_metadata %<>% dplyr::select(-c(PatientID))
  #}
  }
  
  if (variable=="Group"){
    if ("post_LTx" %in% metadata[,variable]) {
      colors <- c("#309f87","#f9c675","#425387")
      metadata[,variable] <- factor(metadata[,variable],levels = c("HC","pre_LTx","post_LTx"))
    }
    else if (("rPSC" %in% metadata[,variable]) & 
             ("pre_LTx" %in% metadata[,variable])) {
      colors <- c("#309f87","#f9c675","#F08080","#A00000")
      metadata[,variable] <- factor(metadata[,variable],levels = c("HC","pre_LTx","non-rPSC","rPSC"))
    }
    else if ("ibd" %in% metadata[,variable]){
      colors <- c("#A06A2C", "#B2182B")  
      metadata[,variable] <- factor(metadata[,variable],levels = c("no_ibd","ibd"))
    } else if ("high" %in% metadata[,variable]){
      colors <- c("#A06A2C", "#B2182B")  
      metadata[,variable] <- factor(metadata[,variable],levels = c("low","high"))
    }
    else {colors <- c("#309f87","#F08080","#A00000")
    metadata[,variable] <- factor(metadata[,variable],levels = c("HC","non-rPSC","rPSC"))
    }
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
  
  normalize_to_range <- function(x) {
    min_x <- min(x)
    max_x <- max(x)
    return(2 * (x - min_x) / (max_x - min_x) - 1)
  }
  
  if (clinical) pca_vec <- apply(pca_vec,2,function(x) normalize_to_range(x))

  x_lab = paste("PCo1 ", "(",round(imp_vec[1]*100,2),"%", ")", sep="")
  y_lab = paste("PCo2 ", "(",round(imp_vec[2]*100,2),"%", ")", sep="")
  
  p <- ggplot(data_for_pca) + 
    geom_point(aes(x=pca_vec[,1],
                   y=pca_vec[,2],
                   color=!!sym(variable),
                   shape=Country),
               show.legend = show_legend,size=size) + 
    xlab(x_lab)+
    ylab(y_lab)+
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    theme(legend.position = "right") 
  
  if (variable=="Group") p <- p + 
    scale_color_manual(values=colors) +
    stat_ellipse(aes(x=pca_vec[,1],y=pca_vec[,2],color=!!sym(variable)),show.legend = FALSE) 
  else p <- p + scale_color_gradient2(low = "lightblue", mid = "blue", high = "black")
  if (clinical){
    set.seed(123)
    if ("PatientID" %in% colnames(clinical_metadata)){
      my_df <- pca_vec  %>% as.data.frame()
      clinical_metadata <- clinical_metadata[rownames(my_df),]
      #my_df <- my_df[rownames(clinical_metadata),]
      my_df <- merge(my_df %>% rownames_to_column("SampleID"),
                     clinical_metadata %>% rownames_to_column("SampleID") %>%
                       dplyr::select(SampleID,PatientID),by="SampleID",all.x=TRUE) %>% 
        as.data.frame() 
      my_df <- my_df %>%
        group_by(PatientID) %>%
        distinct(PatientID, .keep_all = TRUE) %>%
        as.data.frame() %>%
        column_to_rownames("SampleID") %>%
        dplyr::select(-PatientID) 
      
      sample_names <- rownames(my_df)
      ord$vectors <- as.data.frame((lapply(my_df, as.numeric))) %>%
        `rownames<-`(sample_names) %>% as.matrix()
      clinical_metadata %<>% dplyr::select(-PatientID)
      clinical_metadata <- clinical_metadata[rownames(ord$vectors),]
    }
    
    
    fit <- envfit(ord=ord$vectors, env=clinical_metadata, perm = 999,na.rm=TRUE)
    #vector_coordinates <- as.data.frame(scores(fit, "vectors")) #* ordiArrowMul_custom(fit, ord)
    #vector_coordinates <- as.data.frame(fit$vectors$arrows)
    vector_coordinates <- as.data.frame(scores(fit,"vectors")) * ordiArrowMul(fit)
    vector_coordinates <- vector_coordinates[fit$vectors$pvals < 0.05,]
    #factor_coordinates <- as.data.frame(scores(fit, "factors")) #* ordiArrowMul_custom(fit,ord)
    factor_coordinates <- as.data.frame(fit$factors$arrows)
    factor_coordinates <- factor_coordinates[fit$factors$pvals < 0.05,] * ordiArrowMul(fit)
    pca_loadings <- rbind(vector_coordinates,factor_coordinates)
    
    # extracts vector of r2 for each env. variable
    
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


pca_plot_corr <- function(asv_table,taxa_table,metadata, 
                            measure="robust.aitchison",
                            show_boxplots = TRUE,
                            variable = "Group", size=2,
                            segment=NULL,
                            show_legend=TRUE,
                            clinical=FALSE,
                            clinical_metadata=NULL,
                          variables=NULL){
  # Generates a PCA plot showing the beta diversity of the dataset 
  # inputs:
  # asv_table - data frame with column 'SeqID' corresponding to unique ID for ASV or taxa
  # taxa_table - Taxonomy with 'SeqID' and taxonomy
  # metadata - data frame with 'SampleID' identifier
  # measure - beta diversity distance metric, default robust.aitchison
  # show_boxplots - boolean, should boxplots be visualized along the PC, default TRUE 
  # show_legend - boolean, should legend be visualized?, default TRUE 
  # variable - variable that should be used for color-coding, default 'Group'
  # size - size of the points, default 2
  # clinical - boolean, if clinical variables should be visualized as BIPLOT, default FALSE
  # clinical_metadata, if clinical=TRUE, clinical_metadata should be provided with 'SampleID' identifier
  # outputs:
  # p - PCA plot, ggplot object
  
  if (variable=="Group"){
    metadata <- metadata %>%
      mutate(Group = case_when(
        Group == "healthy" ~ "HC",
        Group == "post_ltx" ~ "post_LTx",
        Group == "pre_ltx" ~ "pre_LTx",
        TRUE ~ Group  # Keep other values as is
      ))
  }
  
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
      dplyr::select(-c(Group,Country))
    
    #if (length(unique(clinical_metadata$PatientID))==nrow(clinical_metadata)){
    #  clinical_metadata %<>% dplyr::select(-c(PatientID))
    #}
  }
  
  if ("post_LTx" %in% metadata[,variable]) {
    colors <- c("#309f87","#f9c675","#425387")
    metadata$Group <- factor(metadata$Group,levels = c("HC","pre_LTx","post_LTx"))
  }
  else if (("rPSC" %in% metadata[,variable]) & 
           ("pre_LTx" %in% metadata$Group)) {
    colors <- c("#309f87","#f9c675","#F08080","#A00000")
    metadata$Group <- factor(metadata$Group,levels = c("HC","pre_LTx","non-rPSC","rPSC"))
  }
  else if ("ibd" %in% metadata[,variable]){
    colors <- c("#A06A2C", "#B2182B")  
    metadata$Group <- factor(metadata$Group,levels = c("no_ibd","ibd"))
  } else if ("high" %in% metadata[,variable]){
    colors <- c("#A06A2C", "#B2182B")  
    metadata$Group <- factor(metadata$Group,levels = c("low","high"))
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
  
  normalize_to_range <- function(x) {
    min_x <- min(x)
    max_x <- max(x)
    return(2 * (x - min_x) / (max_x - min_x) - 1)
  }
  
  if (clinical) pca_vec <- apply(pca_vec,2,function(x) normalize_to_range(x))
  
 if (clinical){
    set.seed(123)
    if ("PatientID" %in% colnames(clinical_metadata)){
      my_df <- pca_vec  %>% as.data.frame()
      clinical_metadata <- clinical_metadata[rownames(my_df),]
      my_df <- merge(my_df %>% rownames_to_column("SampleID"),
                     clinical_metadata %>% rownames_to_column("SampleID"),by="SampleID",all.x=TRUE) %>% 
        as.data.frame() 
      my_df <- my_df %>%
        group_by(PatientID) %>%
        distinct(PatientID, .keep_all = TRUE) %>%
        as.data.frame() %>%
        column_to_rownames("SampleID") %>%
        dplyr::select(-PatientID) 
      
      sample_names <- rownames(my_df)
      clinical_metadata %<>% dplyr::select(-PatientID)
      clinical_metadata <- clinical_metadata[rownames(ord$vectors),]
      
  } else {
   my_df <- pca_vec  %>% as.data.frame()
   clinical_metadata <- clinical_metadata[rownames(my_df),]
   my_df <- merge(my_df %>% rownames_to_column("SampleID"),
                  clinical_metadata %>% rownames_to_column("SampleID"),
                  by="SampleID",all.x=TRUE) %>% 
     as.data.frame() 
  }
    
   if (is.null(variables)) variables <- colnames(clinical_metadata)
   
   p_values <- c()
   r_values <- c()
   for (variable in variables){
    cor_res <- cor.test(my_df$Axis.1,my_df[,variable], cor.method="spearman",rm.na=TRUE)
    p_val <- cor_res$p.value
    r_val <- cor_res$estimate
    p_values <- c(p_values,p_val)
    r_values <- c(r_values,r_val)
   }
   p_values_adjusted <- p.adjust(p_values,method="BH")
   names(r_values) <- variables
   r_values[p_values_adjusted<0.05]
   res <- data.frame(variable=variables,
                     r=r_values,
                     p=p_values,
                     p.adj=p_values_adjusted)
   p_values <- c()
   r_values <- c()
   for (variable in variables){
     cor_res <- cor.test(my_df$Axis.2,my_df[,variable], cor.method="spearman",rm.na=TRUE)
     p_val <- cor_res$p.value
     r_val <- cor_res$estimate
     p_values <- c(p_values,p_val)
     r_values <- c(r_values,r_val)
   }
   p_values_adjusted <- p.adjust(p_values,method="BH")
   names(r_values) <- variables
   r_values[p_values_adjusted<0.05]
   res2 <- data.frame(variable=variables,
                     r=r_values,
                     p=p_values,
                     p.adj=p_values_adjusted)
 } 

  if (segment=="terminal_ileum"){
    samples_to_test <- rownames(my_df %>% column_to_rownames("SampleID"))
    mdi_bact <- read.xlsx("../results/Q1/univariate_analysis/psc_effect_terminal_ileum.xlsx")$SeqID
    
    
  } else if (segment=="colon"){
    samples_to_test <- rownames(my_df)
    mdi_bact <- read.xlsx("../results/Q1/univariate_analysis/psc_effect_colon.xlsx")$SeqID
    
  } else (cat("Invalid segment"))
  
  bacteria_df <- vegan::decostand(asv_table %>% column_to_rownames("SeqID"),method = "clr", MARGIN = 2, pseudocount=0.5) %>% as.matrix()
  bacteria_df <- bacteria_df[,samples_to_test] %>% as.data.frame() %>%  rownames_to_column("SeqID")
  bacteria_df %<>% dplyr::filter(SeqID %in% mdi_bact)  %>% column_to_rownames("SeqID")
  pc1 <- my_df$Axis.1
  pc2 <- my_df$Axis.2
  
  corr_coef <- c()
  p_values1 <- c()
  r_values1 <- c()
  p_values2 <- c()
  r_values2 <- c()
  for (bact in rownames(bacteria_df)){
    cor_res <- cor.test(as.vector(as.matrix(bacteria_df[bact,])),pc1,method="spearman")
    r_values1 <- c(r_values1,cor_res$estimate)
    p_values1 <- c(p_values1,cor_res$p.value)
    
    cor_res <- cor.test(as.vector(as.matrix(bacteria_df[bact,])),pc2,method="spearman")
    r_values2 <- c(r_values2,cor_res$estimate)
    p_values2 <- c(p_values2,cor_res$p.value)
  }
  
  names(r_values1) <- rownames(bacteria_df)
  names(p_values1) <- rownames(bacteria_df)
  p_adjusted1 <- p.adjust(p_values1,method = "BH")
  res_mdi1 <- data.frame(r=r_values1,
                         p=p_values1,
                         p_adjusted=p_adjusted1)
  
  names(r_values2) <- rownames(bacteria_df)
  names(p_values2) <- rownames(bacteria_df)
  p_adjusted2 <- p.adjust(p_values2,method = "BH")
  res_mdi2 <- data.frame(r=r_values2,
                         p=p_values2,
                         p_adjusted=p_adjusted2)
  
  #res_mdi1 %<>% filter(p_adjusted<0.05)
  #res_mdi2 %<>% filter(p_adjusted<0.05)
  
  return(list(res,res2,res_mdi1,res_mdi2))
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


heatmap_linda <- function(linda.output,taxa_tab){
  # Generates a heatmap showing the significant taxa in DAA with 
  # the significance information and logfoldchange color-coding
  # inputs:
  # linda.output - output of linda
  # taxa_tab - Taxonomy with 'SeqID', will be used for the row-naming
  # outputs:
  # p - heatmap
  
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
  # Generates a dotheatmap showing the significant taxa in DAA with 
  # the significance information, logfoldchange color-coding and abundance size-coding
  # inputs:
  # raw_linda - univariate statistics for this comparison (uni_statistics variable)
  # taxa_table - Taxonomy with 'SeqID', will be used for the row-naming
  # group - name of the group that should be visualized. If null, all groups will be plotted, default NULL
  # outputs:
  # p - dotheatmap
  
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
  } else if (any(grepl("pre_LTx",plot_df_combined$Variable)) & any(grepl("rPSC",plot_df_combined$Variable))){
    plot_df_combined$Variable <- factor(plot_df_combined$Variable,
                                        levels=c("pre_LTx vs HC","rPSC vs non-rPSC"))
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
  # Generates a volcano showing the significant taxa in linDA 
  # inputs:
  # linda.output - linda output
  # group - name of the comparison that should be plotted
  # taxa_table - Taxonomy with 'SeqID', will be used for the row-naming
  # cutoff.pval - significance threshold, default 0.05
  # cutoff.lfc - logfoldchange threshold, default 1
  # outputs:
  # p - volcano_plot
  
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
  } 
  else if ("Domain" %in% colnames(taxa_table)){
    data_called$name <- taxa_table[data_called$name,"Domain"]
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


volcano_plot_maaslin <- function(maaslin_output,taxa_table,cutoff.pval=0.05, cutoff.lfc=1,variable="Group"){
  # Generates a volcano showing the significant taxa in Maaslin 
  # inputs:
  # maaslin_output - masslin output
  # variable - variable that should be plotted, default "Group"
  # taxa_table - Taxonomy with 'SeqID', will be used for the row-naming
  # cutoff.pval - significance threshold, default 0.05
  # cutoff.lfc - logfoldchange threshold, default 1
  # outputs:
  # p - volcano_plot
  
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


roc_curve <- function(model_object, group){
  # Generates a ROC curve as the performance metric of a model 
  # inputs:
  # model_object - model object
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # outputs:
  # p - roc_curve
  
  ggroc_data <- ggroc(model_object$kfold_rocobjs)$data
  roc_c <- ggplot(data=ggroc_data) + 
    geom_line(aes(x=`1-specificity`, y=sensitivity, by=name, color="red",alpha=0.9)) +
    theme_minimal() + 
    theme(legend.position = "none") + 
    ggtitle(paste0(group[1],' vs ',group[2],
                   ' (AUC = ', round(model_object$model_summary$auc_optimism_corrected,2), 
                   ', P = ',(mean(model_object$valid_performances$auc_validation<0.5)*2),')'))  
    return(roc_c)
}

roc_curve_all <- function(objects){
  # Generates a ROC curve for MDI discriminative power visualization
  # This functions will put all comparisons in one plot
  # inputs:
  # objects - roc_cs objects, where plain ROC information is stored
  # outputs:
  # p - roc_curves plot
  
  p <- ggplot()
  colors <- c(
    "red", "blue","#2B9D2B", "#D90368",
    "#984ea3", "#FAF33E","#DF75BE", "grey",
    "#17BACB", "#66c2a5","#A5BE00", "#000000",
    "#a65629")
  
  colors <- c("#4169E1","#984ea3","#008080",
              "#FF6347","#FFD700","#D90368")
  
  
  names(colors) <- names(objects)
  
  for (i in 1:length(objects)){
    auc <- objects[[i]]$auc
    my_df <- data.frame(Sensitivity=objects[[i]]$sensitivities,
                        `1-specificity`=1-objects[[i]]$specificities,
                        check.names = FALSE)
    
    
    my_color <- names(colors)[i]
    p <- p + 
      geom_line(data=my_df,aes(x=`1-specificity`, y=Sensitivity, group=name,
                               color=!!my_color), linewidth=1.5,alpha=1) 
    
    
    
  }
  p <- p + theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    scale_x_continuous(breaks=c(0,0.5,1)) + 
    scale_y_continuous(breaks=c(0,0.5,1)) + 
    scale_color_manual(values = colors) + #+ guides(colour = "none") + 
    theme(legend.title = element_blank())
  return(p)
}

roc_curve_all_custom <- function(objects,Q,model_name,legend=TRUE){
  # Generates a ROC curve as the performance metric of a model
  # This functions will put all comparisons in one plot
  # inputs:
  # objects - roc_cs objects, where plain ROC information is stored
  # Q - question type, will be used for the path where model is stored (Q1/Q2/Q3)
  # model_name - name of the model
  # legend - boolean, if legend should be generated, default TRUE
  # outputs:
  # p - roc_curves plot
  
  #print(names(objects))
  p <- ggplot()
  
  if (grepl("^Q1", Q)) colors <- c("#FF6347","#708090","#17BACB")
  if (grepl("^Q2",Q)) {
    #colors <- c("#FF7F50","#FFD700","#4169E1","#008080") 
    #colors <- c("#008080","#FFD700","#4169E1","#FF7F50")
    colors <- c("#008080","#FF7F50")
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
    
    my_color <- colors[i]
    
    p <- p + 
      #  geom_line(data=df,aes(x=x,y = y1,color=!!my_color)) +
      geom_ribbon(data=df,aes(x =x,ymin = pmin(y1, y2), ymax = pmax(y1, y2)), fill=my_color, 
                  alpha = 0.5,show.legend = TRUE) +
      #geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1.5) +
      theme_minimal() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            panel.grid = element_blank(),
            axis.ticks.x = element_line(size=0.3,color = "black"),
            axis.ticks.y = element_line(size=0.3,color="black"),
            axis.ticks.length = unit(4,"pt")) + 
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
    
    my_color <- colors[i]
    
    p <- p + 
      geom_line(data=df,aes(x=x, y = y3),color=my_color,size=1.5) +
      theme_minimal() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
            axis.ticks.x = element_line(size=0.3,color = "black"),
            axis.ticks.y = element_line(size=0.3,color="black"),
            axis.ticks.length = unit(4,"pt"), 
            panel.grid = element_blank()) + 
      #scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
      #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
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

## DAA functions ----

group_intersection <- function(group, list_intersections, list_venns,
                               linda.output1, fit_data,
                               raw_linda_results, segment,level){
  # Performs an intersection between linda and maaslin
  # inputs:
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # list_intersections - list for storing the intersection taxa
  # list_venns - list for storing the venn diagrams
  # linda.output1 - linda output
  # fit_data - maaslin output
  # raw_linda_results - linda output with univariate statistics
  # segment - terminal_ileum/colon
  # level - ASV/genus/phylum or other
  # outputs:
  # (list(list_intersections,list_venns,venn))
  
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

country_interaction <- function(group,linda.output1,list_intersections, 
                                uni_data, uni_metadata,
                                segment, level){
  # Finds a significant interaction effect in linDA output and removes them
  # inputs:
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # list_intersections - list for storing the intersection taxa
  # linda.output1 - linda output
  # uni_data - ASV table used in DAA
  # uni_metadata - metadata with 'SampleID' identifier
  # segment - terminal_ileum/colon
  # level - ASV/genus/phylum or other
  # outputs:
  # (list(intersection_significant,linda_czech_interaction_significant,linda_no_interaction_significant))
  # - significant interaction effect in intersection list
  
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
  
  # Removes unsuitable significant taxa from DAA results, only ones with the same direction in 
  # both directions should be retained.
  # inputs:
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # list_intersection_significant - list with problematic taxa, output of country_interaction()
  # list_intersections - original list with the intersection of linda and maaslin, created by group_intersection()
  # segment - terminal_ileum/colon
  # level - ASV/genus/phylum or other
  # outputs:
  # list_intersections - final significant taxa
  
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
rawlinda.df <- function(linda.output,group,uni_data,uni_taxa){
  # creates a dataframe from raw linda output
  # inputs:
  # linda.output - linda output
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # uni_data - dataframe used for DAA
  # uni_taxa - taxonomy with 'SeqID' column
  # outputs:
  # raw_linda_result - dataframe generated from linda.output (all taxa included)
  
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
  # Similar function to rawlinda.df, but here only significant taxa are stored
  # inputs:
  # linda.output - linda output
  # group - names of groups that were compared, e.g. c("rPSC","non-rPSC")
  # uni_data - dataframe used for DAA
  # uni_taxa - taxonomy with 'SeqID' column
  # outputs:
  # linda_result - dataframe generated from linda.output (only significant taxa included)
  
  if (ncol(filt_ileum_uni_taxa)==2){
    seq_ids <- filt_ileum_uni_taxa$SeqID
    filt_ileum_uni_taxa <- data.frame(Domain=filt_ileum_uni_taxa$Domain) %>% `row.names<-`(seq_ids)
  } else filt_ileum_uni_taxa <- filt_ileum_uni_taxa %>% column_to_rownames("SeqID")
  diff_asvs <- rownames(linda.output[[group]])[linda.output[[group]]$padj < 0.05]
  
  # save the results
  if (length(diff_asvs)>0){
    linda_result <- data.frame(
    ASVs=diff_asvs,
    Taxonomy=create_asv_taxa_table(
      filt_ileum_uni_data[diff_asvs,] %>% rownames_to_column("SeqID"),
      filt_ileum_uni_taxa[diff_asvs,] %>% rownames_to_column("SeqID"))$SeqID,
    log2FoldChange=linda.output[[group]][diff_asvs,"log2FoldChange"],
    pvalue=linda.output[[group]][diff_asvs,"pvalue"],
    padj=linda.output[[group]][diff_asvs,"padj"])
  } else linda_result <- NULL
  
  
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

cliffs_delta <- function(data,metadata, groups){
  # This function calculates the cliffs_delta

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
  # Calculates basic univariate statistics (clr as well as ra) for each taxon in the dataset
  # inputs:
  # uni_data - list(asv_table,taxa_table,metadata)
  # outputs:
  # df - dataframe with Q1, median, Q3, prevalence etc.
  
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

  return(res)
}


## Machine learning functions ----
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
  
  # Fits a GLMNET model to pre-prepared dataset by binomial_prep(usage="ml_clr/ra")
  # inputs:
  # data - prepared dataframe using binomial_prep(), contains only two groups
  # outcome - vector of outcome (labels)
  # sample_method - atypboot - out-of-sample boostrap
  # clust_var, name of clustering variable (here, it can be Patient), default NULL
  # N, number of bootstrap samples, default 10 (minimum 100 needed for reportable results)
  # alphas - vector of alpha values to test in glmnet, default seq(0, 1, by = 0.2)
  # family - A string specifying the family for the GLMNET model, default 'binomial'
  # overfitting_check - boolean, if random labels reshuffling should be performed, default FALSE
  # seed - random seed for reproducibility, default 123
  # reuse - boolean, should model just be reloaded?, default FALSE
  # file - name of the file for loading the pre-trained model
  # Q - analysis question - Q1/Q2/Q3..
  # outputs: 
  # enet_model - list(model_summary,valid_performances,
  #                   valid_performance,
  #                   predictions, 
  #                   betas,
  #                   conf_matrices,
  #                   roc_curve,
  #                   kfold_roc_curves,
  #                   calibration_plot, 
  #                   trained_model)
  
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
  # Fits a Gboost model to dataset prepared by binomial_prep(usage="ml_clr/ra")
  # inputs:
  # data - prepared dataframe using binomial_prep(), contains only two groups
  # outcome - vector of outcome (labels)
  # sample_method - atypboot - out-of-sample boostrap
  # clust_var, name of clustering variable (here, it can be Patient), default NULL
  # N, number of bootstrap samples, default 10 (minimum 100 needed for reportable results)
  # family - A string specifying the family for the GLMNET model, default 'binomial'
  # overfitting_check - boolean, if random labels reshuffling should be performed, default FALSE
  # seed - random seed for reproducibility, default 123
  # reuse - boolean, should model just be reloaded?, default FALSE
  # file - name of the file for loading the pre-trained model
  # Q - analysis question - Q1/Q2/Q3..
  # outputs: 
  # gbm_model - list(model_summary,valid_performances,
  #                   valid_performance,
  #                   predictions, 
  #                   roc_curve,
  #                   kfold_roc_curves,
  #                   trained_model)
  
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
  
  # Fits a kNN model to dataset prepared by binomial_prep(usage="ml_clr/ra")
  # inputs:
  # data - prepared dataframe using binomial_prep(), contains only two groups
  # outcome - vector of outcome (labels)
  # sample_method - atypboot - out-of-sample boostrap
  # clust_var, name of clustering variable (here, it can be Patient), default NULL
  # N, number of bootstrap samples, default 10 (minimum 100 needed for reportable results)
  # family - A string specifying the family for the GLMNET model, default 'binomial'
  # overfitting_check - boolean, if random labels reshuffling should be performed, default FALSE
  # seed - random seed for reproducibility, default 123
  # reuse - boolean, should model just be reloaded?, default FALSE
  # file - name of the file for loading the pre-trained model
  # Q - analysis question - Q1/Q2/Q3..
  # outputs: 
  # knn_model - list(model_summary,valid_performances,
  #                   valid_performance,
  #                   predictions, 
  #                   roc_curve,
  #                   kfold_roc_curves,
  #                   trained_model)
  
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
                         family = 'binomial',
                         overfitting_check=FALSE,
                         seed = 123,
                         reuse=FALSE,
                         file=NULL,
                         Q=NULL) {

  # Fits a Gboost model to dataset prepared by binomial_prep(usage="ml_clr/ra")
  # inputs:
  # data - prepared dataframe using binomial_prep(), contains only two groups
  # outcome - vector of outcome (labels)
  # sample_method - atypboot - out-of-sample boostrap
  # clust_var, name of clustering variable (here, it can be Patient), default NULL
  # N, number of bootstrap samples, default 10 (minimum 100 needed for reportable results)
  # family - A string specifying the family for the GLMNET model, default 'binomial'
  # overfitting_check - boolean, if random labels reshuffling should be performed, default FALSE
  # seed - random seed for reproducibility, default 123
  # reuse - boolean, should model just be reloaded?, default FALSE
  # file - name of the file for loading the pre-trained model
  # Q - analysis question - Q1/Q2/Q3..
  # outputs: 
  # rf_model - list(model_summary,valid_performances,
  #                   valid_performance,
  #                   predictions, 
  #                   roc_curve,
  #                   kfold_roc_curves,
  #                   trained_model)
  
  if (all(data[,1] >= 0 & data[,1] <= 1)) ra = TRUE
  else  ra = FALSE
  
  if (overfitting_check) {
    data$Group <- sample(data$Group)
  }
  
  if (reuse){
    if (ra) {
      if (overfitting_check) load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model_ra.RData"))
      else load(file.path("../intermediate_files/models/",Q,file,"rf_model_ra.RData"))
    } else {
      if (overfitting_check) {
        load(file.path("../intermediate_files/models_overfitting_check/",Q,file,"rf_model.RData"))
      }
      else load(file.path("../intermediate_files/models/",Q,file,"rf_model.RData"))
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

# Clinical analysis ------------------

ordiArrowMul_custom <- function (x, ord, at = c(0,0), fill = 0.75,
                                 display, choices = c(1,2)) {
  if (length(x$factors)>0) X <- do.call(rbind,scores(x,c("vectors", "factors")))
  else X <- scores(x,c("vectors"))
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


dysbiosis_index_calculation <- function(my_table, metadata_table, 
                                        psc_increased,
                                        psc_decreased,
                                        name){
  # Calculates the microbial dysbiosis index (MDI) on relative abundances
  # inputs:
  # my_table - ASV table with 'SeqID'
  # metadata_table - metadadata with 'SampleID'
  # psc_increased - names of increased taxa in PSC
  # psc_descread - names of decreased taxa in PSCA
  # name - string indicating which index to calculate (e.g. dys_unfiltered_asv)
  # output: 
  # dysbiosis_data - dataframe with MDI column
  
  my_table <- my_table %>% column_to_rownames("SeqID")
  
  dysbiosis_data <- data.frame()
  for (i in 1:ncol(my_table)){
    SampleID <- colnames(my_table)[i]
    PatientID <- metadata_table[metadata_table$SampleID==SampleID,"Patient"]
    abundances <- my_table[,i]/sum(my_table[,i])
    names(abundances) <- rownames(my_table)
    abundances_psc_increased <- sum(abundances[psc_increased])
    abundances_psc_decreased <- sum(abundances[psc_decreased])
    
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
  # Similar function to dysbiosis_index_calculation(), but this is for CLR-trasformed data
  
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


clinical_boxplot <- function(metadata,variable){
  # Creates boxplot for clinical variables
  # inputs:
  # metadata - metadata dataframe with 'SampleID'
  # variable - name of variable to be plotted
  
  if (!("PSC_IBD" %in% colnames(metadata))){
    metadata_var <- metadata %>% 
      dplyr::select(all_of(c("PatientID",variable, "Group","Country"))) %>% 
      as.data.frame()
    
    metadata_var[,variable] <- as.numeric(metadata_var[,variable])
    metadata_var_without_na <- metadata_var[!is.na(metadata_var[,variable]),]
    
    colors <- c("#309f87","#f9c675","#F08080","#A00000")
    if (length(unique(metadata_var_without_na$Group))==1) colors <- c("#f9c675")
    if (!"healthy" %in% tolower(unique(metadata_var_without_na$Group))) {
      if ("rpsc" %in% tolower(unique(metadata_var_without_na$Group))) colors <- c("#f9c675","#F08080","#A00000")
      else colors <- c("#f9c675","#425387")
    }  else if (("rpsc" %in% tolower(unique(metadata_var_without_na$Group))) &
             (!"pre_ltx" %in% tolower(unique(metadata_var_without_na$Group)))) {
      colors <- c("#309f87","#F08080","#A00000")
    } else if (("pre_ltx" %in% tolower(unique(metadata_var_without_na$Group))) &
               ("post_ltx" %in% tolower(unique(metadata_var_without_na$Group)))) {
      colors <- c("#309f87","#f9c675","#425387")
    }
  } else {
    metadata_var <- metadata %>% dplyr::select(all_of(c("PatientID",
                                                       variable, "PSC_IBD","Country"))) %>% 
      as.data.frame() %>%
      dplyr::mutate(Group=PSC_IBD)
    metadata_var[,variable] <- as.numeric(metadata_var[,variable])
    metadata_var_without_na <- metadata_var[!is.na(metadata_var[,variable]),]
    colors <- c("#A06A2C", "#B2182B") 
  }
  
  
  p <- ggplot(metadata_var_without_na) + 
    geom_boxplot(aes(x=Group, y=!!sym(variable)),outliers = FALSE) + 
    geom_jitter(width = 0.2,height = 0,aes(x=Group, y=!!sym(variable), color=Group,shape=Country),size=2) +
    scale_fill_manual(values=colors) + 
    scale_color_manual(values=colors) + 
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank())
  return(p)
}



clinical_correlation <- function(metadata,variable,level,
                                 segment="ileum"){
  # Calculates spearman correlations between clinical variables and MDI
  # inputs:
  # metadata - metadata dataframe with 'SampleID' and 'MDI'
  # variable - name of clinical variable - column in metadata df
  # level - ASV/genus, on which level MDI was calculated
  # segment - terminal_ileum, colon - in colon, 100 iterations will be performed
  # outputs:
  # cor - table with $P (p_values) and $r (correlation coefficients)
  
  set.seed(123)
  level <- tolower(level)
  
  metadata_var <- metadata %>% 
    dplyr::select(all_of(c("SampleID", "PatientID",
                           paste0("dys_filtered_",level),
                           variable, "Group")))
  
  metadata_var[,variable] <- as.numeric(metadata_var[,variable])
  metadata_var <- drop_na(metadata_var)
  
  if (segment=="colon"){
    n_iterations <- 100
    corrs <- list()
    for (i in 1:n_iterations){
      random_sampled_df <- metadata_var %>%
        group_by(PatientID) %>%  # Group by PatientID
        slice_sample(n = 1) %>% # Randomly pick one row per PatientID
        ungroup()
      
      corr <- cor.table(
        random_sampled_df[,c(variable,paste0("dys_filtered_",level))], 
        cor.method="spearman"
      )
      
      corr$P <- corr$P[2]
      corr$r <- corr$r[2]
      corrs[[i]] <- corr
    }
    
    r_values_list <- unlist(map(corrs, ~ .x$r))  # List of r matrices
    p_values_list <- unlist(map(corrs, ~ .x$P))  # List of P matrices
    mean_r <- round(median(r_values_list),2)
    
    p_value <- quantile(p_values_list,0.9)
    corr <- list()
    corr$r <- mean_r
    corr$P <- p_value
  } else {
    corr <- cor.table(
      metadata_var[,c(variable,paste0("dys_filtered_",level))],
      cor.method="spearman")
    
    corr$P <- corr$P[2]
    corr$r <- round(corr$r[2],2)
  }
  
  return(corr)
}

clinical_correlation_abundances <- function(metadata,variable,taxon,level,
                                            segment="ileum",rename_p=FALSE){
  # Calculates spearman correlations between taxon and MDI
  # inputs:
  # metadata - metadata dataframe with 'SampleID' and 'MDI'
  # variable - name of clinical variable - column in metadata df
  # taxon - name of taxon
  # level - ASV/genus, on which level MDI was calculated
  # segment - terminal_ileum, colon - in colon, 100 iterations will be performed
  # outputs:
  # cor - table with $P (p_values) and $r (correlation coefficients)
  
  set.seed(123)
  level <- tolower(level)
  
  metadata_var <- metadata %>% dplyr::select(all_of(c("SampleID", "PatientID",
                                                      variable,taxon,
                                                      "Group")))
  metadata_var[,variable] <- as.numeric(metadata_var[,variable])
  metadata_var <- drop_na(metadata_var)
  
  if (segment=="colon"){
    n_iterations <- 100
    corrs <- list()
    for (i in 1:n_iterations){
      random_sampled_df <- metadata_var %>%
        group_by(PatientID) %>%  # Group by PatientID
        slice_sample(n = 1) %>% # Randomly pick one row per PatientID
        ungroup()
      
      corr <- cor.table(random_sampled_df[,c(variable,taxon)], cor.method="spearman")
      corr$P <- corr$P[2]
      corr$r <- corr$r[2]
      corrs[[i]] <- corr
    }
    r_values_list <- unlist(map(corrs, ~ .x$r))  # List of r matrices
    p_values_list <- unlist(map(corrs, ~ .x$P))  # List of P matrices
    mean_r <- round(median(r_values_list),2)
    
    p_value <- quantile(p_values_list,0.9)
    if (rename_p) p_value <- ifelse(p_value < 0.001, " < 0.001", ifelse(p_value < 0.01, " < 0.01", ifelse(p_value < 0.05, " < 0.05", paste0("=",round(p_value,3)))))
    corr <- list()
    corr$r <- mean_r
    corr$P <- p_value
  } else {
    corr <- cor.table(metadata_var[,c(variable,taxon)], cor.method="spearman")
    corr$r <- round(corr$r[2],2)
    p_value <- corr$P[2]
    if (rename_p) p_value <- ifelse(p_value < 0.001, " < 0.001", ifelse(p_value < 0.01, " < 0.01", ifelse(p_value < 0.05, " < 0.05", paste0("=",round(p_value,3)))))
    corr$P <- p_value
  }
  
  return(corr)
}

clinical_scatter <- function(corr,metadata,variable,level){
  # Creates scatterplot visualizing MDI~clinical_variable and its correlation 
  # coefficient
  # inputs:
  # corr - correlations, calculated by clinical_correlation()
  # metadata - metadata dataframe with 'SampleID'
  # variable - name of variable to be plotted
  # level - ASV/genus ...
  
  level <- tolower(level)
  df_for_scatter_plot <- metadata %>% 
    dplyr::select(all_of(c("SampleID", "PatientID",
                           paste0("dys_filtered_",level),
                           variable, "Group","Country")))
  
  df_for_scatter_plot[,variable] <- as.numeric(df_for_scatter_plot[,variable])
  
  df_for_scatter_plot <- drop_na(df_for_scatter_plot)
  
  colors <- c("#309f87","#f9c675","#F08080","#A00000")
  if (!"healthy" %in% unique(df_for_scatter_plot$Group)) colors <- c("#f9c675","#F08080","#A00000")
  
  r_corr <- corr$r
  p_corr <- corr$P
  p_corr <- ifelse(p_corr < 0.001, " < 0.001", ifelse(p_corr < 0.01, " < 0.01", ifelse(p_corr < 0.05, " < 0.05", paste0("=",round(p_corr,3)))))
  
  x_var <- paste0("dys_filtered_",level)
  x_position <- (max(df_for_scatter_plot[,x_var])) - 0.2*(max(df_for_scatter_plot[,x_var]) - min(df_for_scatter_plot[,x_var]))
  
  y_position <- max(df_for_scatter_plot[,variable]) - 0.05*(max(df_for_scatter_plot[,variable]) - min(df_for_scatter_plot[,variable]))
  
  labels <- paste0("R = ",r_corr, "\\nP ",p_corr)
  p <- ggplot(df_for_scatter_plot) + 
    geom_point(aes(x=!!sym(x_var), y=!!sym(variable), color=Group,shape=Country)) + 
    geom_smooth(aes(x=!!sym(x_var), y=!!sym(variable)),method=lm, se=TRUE) +
    geom_text(aes(x=x_position, y=y_position),
              label=gsub("\\\\n", "\n", labels),hjust = 0, size = 5/.pt) + 
    scale_color_manual(values=colors) + 
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) + 
    xlab("MDI")
  return(p)
}

clinical_scatter_abundances <- function(corr,metadata,variable,taxon,level,
                                        size=4){
  # Creates scatterplot visualizing variable~taxon and its correlation 
  # coefficient
  # inputs:
  # corr - correlations, calculated by clinical_correlation()
  # metadata - metadata dataframe with 'SampleID'
  # variable - name of variable to be plotted
  # taxon - name of the taxon to be plotted
  # level - ASV/genus ...
  # size - size of point in ggplot
  # outputs:
  # p - scatterplot 
  
  level <- tolower(level)
  df_for_scatter_plot <- metadata %>% dplyr::select(all_of(c("SampleID", "PatientID","Country",
                                                             taxon,
                                                             variable, "Group")))
  df_for_scatter_plot[,variable] <- as.numeric(df_for_scatter_plot[,variable])
  df_for_scatter_plot <- drop_na(df_for_scatter_plot)
  colors <- c("#309f87","#f9c675","#F08080","#A00000")
  if (!"healthy" %in% unique(df_for_scatter_plot$Group)) colors <- c("#f9c675","#F08080","#A00000")
  
  x_var <- taxon
  r_corr <- corr$r
  p_corr <- corr$P
  
  x_position <- (max(df_for_scatter_plot[,x_var])) - 0.2*(max(df_for_scatter_plot[,x_var]) - min(df_for_scatter_plot[,x_var]))
  y_position <- max(df_for_scatter_plot[,variable]) - 0.05*(max(df_for_scatter_plot[,variable]) - min(df_for_scatter_plot[,variable]))
  labels <- paste0("r = ",r_corr, "\\np ",p_corr)
  p <- ggplot(df_for_scatter_plot) + 
    geom_point(aes(x=!!sym(x_var), y=!!sym(variable), color=Group,shape=Country)) + 
    geom_smooth(aes(x=!!sym(x_var), y=!!sym(variable)),method=lm, se=TRUE) +
    geom_text(aes(x=x_position, y=y_position,label=gsub("\\\\n", "\n", labels)),size=size) + 
    scale_color_manual(values=colors) + 
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0),
          axis.ticks.x = element_line(size=0.3,color = "black"),
          axis.ticks.y = element_line(size=0.3,color="black"),
          axis.ticks.length = unit(4,"pt"),
          panel.grid = element_blank()) 
  return(p)
}



subprepare_for_heatmap <- function(corrs_segment,MDI=TRUE){
  # Function for heatmap construction
  # input:
  # corrs_segment - correlation coefficients calculated for specific segment
  # MDI - boolean, if heatmap for MDI will be created, default TRUE
  # output:
  # list(p_df_sig,r_df) - dataframe with significant correlations to be visualized
  
  corrs_segment <- corrs_segment[!grepl("_log_",names(corrs_segment))]
  names(corrs_segment) <- gsub("score_","",names(corrs_segment))
  # Example list structure for demonstration
  # Initialize empty lists for r and P values
  r_values <- list()
  p_values <- list()
  
  # Extract segment and clinical variables
  for (name in names(corrs_segment)) {
    if (MDI) {
      clinical_var <- sub("^(.*?_){2}", "",  name)
      if (all(grepl("ileum",names(corrs_segment)))) {
        segment <- "MDI terminal ileum"
      } else if (all(grepl("colon",names(corrs_segment)))) {
        segment <- "MDI colon"
      } else (message("PROBLEM with segment"))
      
    }
    else {
      segment <- sub("^(.*?_){3}", "",  name)
      clinical_var <- sub(segment,"", sub("^(.*?_){2}", "",  name))
      clinical_var <- gsub("_","",clinical_var)
    }
    
    
    # Fill the r and P lists
    if (!is.null(corrs_segment[[name]]$r)) {
      r_values[[segment]][clinical_var] <- corrs_segment[[name]]$r
    }
    
    if (!is.null(corrs_segment[[name]]$P)) {
      p_values[[segment]][clinical_var] <- corrs_segment[[name]]$P
    }
  }
  
  # Convert the lists to data frames
  r_df <- do.call(rbind, lapply(r_values, function(x) setNames(as.data.frame(t(x)), names(x))))
  p_df <- do.call(rbind, lapply(p_values, function(x) setNames(as.data.frame(t(x)), names(x))))
  
  # Add row names as segment names
  rownames(r_df) <- names(r_values)
  rownames(p_df) <- names(p_values)
  
  if (MDI) pd_df_corrected <- as.data.frame(p.adjust(p_df,method="BH"))
  else pd_df_corrected <- as.data.frame(apply(p_df,2,function(x) p.adjust(x, method="BH")))
  p_df_sig <- pd_df_corrected
  p_df_sig[,] <- ""
  
  p_df_sig[pd_df_corrected < 0.05] <- "*"
  p_df_sig[pd_df_corrected < 0.01] <- "**"
  p_df_sig[pd_df_corrected < 0.001] <- "***"
  
  if (MDI) p_df_sig <- t(p_df_sig) %>% 
    as.data.frame()  %>% 
    `rownames<-`(segment)
  
  return(list(p_df_sig,r_df))
}

prepare_for_heatmap <- function(corrs_ileum,corrs_colon){
  # Function for heatmap construction
  # input:
  # corrs_ileum - correlation coefficients calculated for terminal_ileum
  # corrs_colon - correlation coefficients calculated for colon
  # output:
  # list(p_df_sig,r_df) - dataframe with significant correlations to be visualized
  
  ileum_list <- subprepare_for_heatmap(corrs_ileum)
  colon_list <- subprepare_for_heatmap(corrs_colon)
  
  p_df_sig <- rbind(ileum_list[[1]],colon_list[[1]])
  r_df <- rbind(ileum_list[[2]],colon_list[[2]])
  
  return(list(p_df_sig,r_df))
  
}
