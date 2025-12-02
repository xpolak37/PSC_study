setwd("~/IKEM/Projects/PSC_study/analysis/results")

compute_se <- function(ci) {
  return((ci[2] - ci[1]) / (2 * 1.96))
}

compare_auc <- function(mean1, se1, mean2, se2) {
  z <- (mean1 - mean2) / sqrt(se1^2 + se2^2)
  p_value <- 2 * (1 - pnorm(abs(z)))  # Two-tailed test
  return(c(z, p_value))
}

run_z_test <- function(groups,models,cohort1,cohort2){
  p_values <- c()
  for (group in groups){
    for (model in models){
      for (segment in c("TI","col")){
        a <- df[paste(group,cohort1),paste(segment,model,sep=".")]
        b <- df[paste(group,cohort2),paste(segment,model,sep=".")]
        
        a <- unlist(strsplit(a,";|\\(|\\)"))
        b <- unlist(strsplit(b,";|\\(|\\)"))
        mean_a <- as.numeric(a[1])
        mean_b <- as.numeric(b[1])
        
        se_a <- compute_se(ci = as.numeric(a[2:3]))
        se_b <- compute_se(ci = as.numeric(b[2:3]))
        
        z_p <- compare_auc(mean_a, se_a, mean_b, se_b)
        p_values <- c(p_values,z_p[2])
        names(p_values)[length(p_values)] <- paste(group,model,segment,sep="_")
      }
    }
  }
  
  p_values_ajusted <- p.adjust(p_values,method = "BH")
  return(list(p_values_ajusted,p_values))
}

run_roc_test <- function(groups,models,cohort1,cohort2){
  p_values <- c()
  std_diffs <- c()
  mean_diffs <- c()
  plot_list <- list()
  for (group in groups){
    group <- gsub("HC","healthy",group)
    for (model in c("knn_model","rf_model","gbm_model","enet_model")){
      for (segment in c("TI","col")){
        if (segment=="TI") file <- "genus terminal_ileum"
        else file <- "genus colon"
        if (grepl("PSC",group)) {
          if (cohort1 == "NOR+CZ") Q_1 <- "Q2"
          else if (cohort1 == "NOR") Q_1 <- "Q2_norway"
          else if (cohort1 == "CZ") Q_1 <- "Q2_czech"
          
          if (cohort2 == "NOR+CZ") Q_2 <- "Q2"
          else if (cohort2 == "NOR") Q_2 <- "Q2_norway"
          else if (cohort2 == "CZ") Q_2 <- "Q2_czech"
          
        } else {
          if (cohort1 == "NOR+CZ") Q_1 <- "Q1"
          else if (cohort1 == "NOR") Q_1 <- "Q1_norway"
          else if (cohort1 == "CZ") Q_1 <- "Q1_czech"
          
          if (cohort2 == "NOR+CZ") Q_2 <- "Q1"
          else if (cohort2 == "NOR") Q_2 <- "Q1_norway"
          else if (cohort2 == "CZ") Q_2 <- "Q1_czech"
        }
        
        # MODELS differences
        model_cohort1 <- get(load(paste0("../intermediate_files/models/",Q_1,"/",group," ", file, paste0("/",model,".RData"))))
        model_cohort2 <- get(load(paste0("../intermediate_files/models/",Q_2,"/",group," ", file, paste0("/",model,".RData"))))
        
        roc1 <- model_cohort1$kfold_rocobjs
        roc2 <- model_cohort2$kfold_rocobjs
      
        
        aucs1 <- sapply(roc1,function(x) x$auc) 
        aucs2 <- sapply(roc2,function(x) x$auc) 
        
        # STD
        # Calculate the differences between model1 and model2
        differences <- aucs1 - aucs2
        
        # Calculate the standard deviation of these differences
        std_diff <- sd(differences)
        std_diffs <- c(std_diffs,std_diff)
        names(std_diffs)[length(std_diffs)] <- paste(group,model,segment,sep="_")
        
        mean_diff <- mean(abs(differences))
        mean_diffs <- c(mean_diffs,mean_diff)
        names(mean_diffs)[length(mean_diffs)] <- paste(group,model,segment,sep="_")
        
        #wilcoxon_result <- wilcox.test(aucs1, aucs2, paired = FALSE)
        wilcoxon_result <- t.test(aucs1, aucs2, paired = FALSE)
        p_values <- c(p_values,wilcoxon_result$p.value)
        names(p_values)[length(p_values)] <- paste(group,model,segment,sep="_")  
        plot_list[[length(plot_list)+1]] <- ggplot(melt(data.frame(`a`=aucs1,`b`=aucs2))) + geom_boxplot(aes(x=variable,y=value))
        }
      }
  }
  p_values_adjusted <- p.adjust(p_values,method="BH")
  return(list(p_values_adjusted, plot_list,std_diffs,mean_diffs))
}

run_z_test_rpsc <- function(groups,models){
  p_values <- c()
    for (model in models){
      for (segment in c("TI","col")){
        a <- df[paste(groups[1],"NOR+CZ"),paste(segment,model,sep=".")]
        b <- df[paste(groups[2],"NOR+CZ"),paste(segment,model,sep=".")]
        
        a <- unlist(strsplit(a,";|\\(|\\)"))
        b <- unlist(strsplit(b,";|\\(|\\)"))
        mean_a <- as.numeric(a[1])
        mean_b <- as.numeric(b[1])
        
        se_a <- compute_se(ci = as.numeric(a[2:3]))
        se_b <- compute_se(ci = as.numeric(b[2:3]))
        
        z_p <- compare_auc(mean_a, se_a, mean_b, se_b)
        p_values <- c(p_values,z_p[2])
        names(p_values)[length(p_values)] <- paste("rPSC eff",model,segment,sep="_")
      }
  }
  
  p_values_ajusted <- p.adjust(p_values,method = "BH")
  return(list(p_values_ajusted,p_values))
}

run_z_test_models <- function(){
  p_values <- c()
  for (group in groups){
    for (model in models){
      for (segment in c("TI","col")){
        a <- df[paste(group,cohort1),paste(segment,model,sep=".")]
        b <- df[paste(group,cohort2),paste(segment,model,sep=".")]
        
        a <- unlist(strsplit(a,";|\\(|\\)"))
        b <- unlist(strsplit(b,";|\\(|\\)"))
        mean_a <- as.numeric(a[1])
        mean_b <- as.numeric(b[1])
        
        se_a <- compute_se(ci = as.numeric(a[2:3]))
        se_b <- compute_se(ci = as.numeric(b[2:3]))
        
        z_p <- compare_auc(mean_a, se_a, mean_b, se_b)
        p_values <- c(p_values,z_p[2])
        names(p_values)[length(p_values)] <- paste(group,model,segment,sep="_")
      }
    }
  }
  
  p_values_ajusted <- p.adjust(p_values,method = "BH")
  return(list(p_values_ajusted,p_values))
  
}

library(openxlsx)
library(pROC)
library(PMCMRplus)
library(reshape2)
library(ggpubr)

df <- read.xlsx("models_all.xlsx",rowNames = TRUE)

groups <- unique(gsub("(NOR)|(CZ)|\\+","",rownames(df)))
groups <- gsub(" $","",groups)
groups1 <- groups[1:3]
groups2 <- groups[4:7]
groups3 <- groups[6:7]

models <- c("knn","rf","gbm","enet")

# countries differences - Q1
## NOR vs CZ
p_nor_cz <- run_z_test(groups1,models,"NOR","CZ")
any(p_nor_cz[[1]] < 0.05,na.rm = TRUE)

p_nor_cz <- run_roc_test(groups1,models,"NOR","CZ")
p_nor_cz
std <- p_nor_cz[[3]]
means <- p_nor_cz[[4]]


# pre_LTx vs HEALTHY TI
mean(std[grepl("(healthy)",names(p_nor_cz[[3]])) & 
    grepl("(pre_LTx)",names(p_nor_cz[[3]])) & 
    grepl("(_TI)$",names(p_nor_cz[[3]])) ])

mean(std[grepl("(healthy)",names(p_nor_cz[[3]])) & 
           grepl("(post_LTx)",names(p_nor_cz[[3]])) & 
           grepl("(_TI)$",names(p_nor_cz[[3]])) ])

mean(means[grepl("(healthy)",names(p_nor_cz[[4]])) & 
           grepl("(pre_LTx)",names(p_nor_cz[[4]])) & 
           grepl("(_TI)$",names(p_nor_cz[[4]])) ])

mean(means[grepl("(healthy)",names(p_nor_cz[[4]])) & 
           grepl("(post_LTx)",names(p_nor_cz[[4]])) & 
           grepl("(_TI)$",names(p_nor_cz[[4]])) ])


# pre_LTx vs HEALTHY COLON
mean(std[grepl("(healthy)",names(p_nor_cz[[3]])) & 
           grepl("(pre_LTx)",names(p_nor_cz[[3]])) & 
           grepl("(_col)$",names(p_nor_cz[[3]])) ])

mean(std[grepl("(healthy)",names(p_nor_cz[[3]])) & 
           grepl("(post_LTx)",names(p_nor_cz[[3]])) & 
           grepl("(_col)$",names(p_nor_cz[[3]])) ])


mean(means[grepl("(healthy)",names(p_nor_cz[[4]])) & 
             grepl("(pre_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(_col)$",names(p_nor_cz[[4]])) ])

mean(means[grepl("(healthy)",names(p_nor_cz[[4]])) & 
             grepl("(post_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(_col)$",names(p_nor_cz[[4]])) ])

# PRE VS POST ILEUM
mean(means[grepl("(post_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(pre_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(_TI)$",names(p_nor_cz[[4]])) ])

# PRE VS POST COLON
mean(means[grepl("(post_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(pre_LTx)",names(p_nor_cz[[4]])) & 
             grepl("(_col)$",names(p_nor_cz[[4]])) ])

## NOR vs NOR+CZ
p_nor_norcz <- run_z_test(groups1,models,"NOR","NOR+CZ")
any(p_nor_norcz[[1]] < 0.05,na.rm = TRUE)

p_nor_norcz <- run_roc_test(groups1,models,"NOR","NOR+CZ")
p_nor_norcz

## CZ vs NOR+CZ
p_cz_norcz <- run_z_test(groups1,models,"CZ","NOR+CZ")
any(p_cz_norcz[[1]] < 0.05,na.rm = TRUE)

p_cz_norcz <- run_roc_test(groups1,models,"CZ","NOR+CZ")
p_cz_norcz

# countries differences - Q2
## NOR vs CZ
p_nor_cz <- run_z_test(groups2,models,"NOR","CZ")
any(p_nor_cz[[1]] < 0.05,na.rm = TRUE)

p_nor_cz <- run_roc_test(groups2,models,"NOR","CZ")
p_nor_cz

## NOR vs NOR+CZ
p_nor_norcz <- run_z_test(groups2,models,"NOR","NOR+CZ")
any(p_nor_norcz[[1]] < 0.05,na.rm = TRUE)

## CZ vs NOR+CZ
p_cz_norcz <- run_z_test(groups2,models,"CZ","NOR+CZ")
any(p_cz_norcz[[1]] < 0.05,na.rm = TRUE)

# rPSC effect differences
p_rPSC <- run_z_test_rpsc(groups3,models)
any(p_rPSC[[1]] < 0.05,na.rm = TRUE)
p_rPSC


for (group in groups){
    group <- gsub("HC","healthy",group)
      for (segment in c("TI","col")){
        if (segment=="TI") file <- "genus terminal_ileum"
        else file <- "genus colon"
        
        if (grepl("PSC",group)) Q <- "Q2"
        else Q <- "Q1"
        # MODELS differences
        load(paste0("../intermediate_files/models/",Q,"/",group," ", file, "/knn_model.RData"))
        load(paste0("../intermediate_files/models/",Q,"/",group," ", file, "/rf_model.RData"))
        load(paste0("../intermediate_files/models/",Q,"/",group," ", file, "/gbm_model.RData"))
        load(paste0("../intermediate_files/models/",Q,"/",group," ", file, "/enet_model.RData"))
        
        roc1 <- knn_model$kfold_rocobjs
        roc2 <- rf_model$kfold_rocobjs
        roc3 <- gbm_model$kfold_rocobjs
        roc4 <- enet_model$kfold_rocobjs
        
        aucs1 <- sapply(roc1,function(x) x$auc)
        aucs2 <- sapply(roc2,function(x) x$auc)
        aucs3 <- sapply(roc3,function(x) x$auc)
        aucs4 <- sapply(roc4,function(x) x$auc)
        
        auc_data <- matrix(data=c(aucs1,aucs2,aucs3,aucs4),ncol=4)
        melted_auc_data <- melt(auc_data)
        colnames(melted_auc_data) <- c("Bootstrap_Iteration","Model","AUC")
        colnames(auc_data) <- models
        res <- friedman.test(AUC ~ Model | Bootstrap_Iteration, data = melted_auc_data)
        
        if (res$p.value < 0.05) {
          print(paste(group,segment))
          posthoc_res <- PMCMRplus::frdAllPairsConoverTest(melted_auc_data$AUC, 
                                                           groups=melted_auc_data$Model,
                                                           blocks=melted_auc_data$Bootstrap_Iteration,
                                                           p.adjust.method = "BH")
          print(posthoc_res)
          
        }
        
      }
  }
        
df_cleaned <- df
df_cleaned[] <- lapply(df, function(x) gsub("\\s?\\(.*\\)", "", x))
df_cleaned_ti <- df_cleaned[,1:4]
df_cleaned_col <- df_cleaned[,5:8]
apply(df_cleaned_ti,1,which.min)
apply(df_cleaned_col,1,which.min)


