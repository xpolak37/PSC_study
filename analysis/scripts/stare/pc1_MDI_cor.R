# colon
samples_to_test <- rownames(my_df)
mdi_bact <- read.xlsx("../results/Q1/univariate_analysis/psc_effect_colon.xlsx")$SeqID
bacteria_df <- vegan::decostand(filt_colon_genus_tab %>% column_to_rownames("SeqID"),method = "clr", MARGIN = 2, pseudocount=0.5) %>% as.matrix()
bacteria_df <- bacteria_df[,rownames(my_df)] %>% as.data.frame() %>%  rownames_to_column("SeqID")
bacteria_df %<>% dplyr::filter(SeqID %in% mdi_bact)  %>% column_to_rownames("SeqID")
pc1 <- my_df$Axis.1

corr_coef <- c()
p_values <- c()
r_values <- c()
for (bact in rownames(bacteria_df)){
  cor_res <- cor.test(as.vector(as.matrix(bacteria_df[bact,])),pc1,method="spearman")
  r_values <- c(r_values,cor_res$estimate)
  p_values <- c(p_values,cor_res$p.value)
}

names(r_values) <- rownames(bacteria_df)
names(p_values) <- rownames(bacteria_df)
p_adjusted <- p.adjust(p_values,method = "BH")
res <- data.frame(r=r_values,
                  p=p_values,
                  p_adjusted=p_adjusted)

res %>% filter(p_adjusted<0.05)

# TI
samples_to_test <- rownames(my_df %>% column_to_rownames("SampleID"))
mdi_bact <- read.xlsx("../results/Q1/univariate_analysis/psc_effect_terminal_ileum.xlsx")$SeqID
bacteria_df <- vegan::decostand(filt_ileum_genus_tab %>% column_to_rownames("SeqID"),method = "clr", MARGIN = 2, pseudocount=0.5) %>% as.matrix()
bacteria_df <- bacteria_df[,samples_to_test] %>% as.data.frame() %>%  rownames_to_column("SeqID")
bacteria_df %<>% dplyr::filter(SeqID %in% mdi_bact)  %>% column_to_rownames("SeqID")
pc1 <- my_df$Axis.2

corr_coef <- c()
p_values <- c()
r_values <- c()
for (bact in rownames(bacteria_df)){
  cor_res <- cor.test(as.vector(as.matrix(bacteria_df[bact,])),pc1,method="spearman")
  r_values <- c(r_values,cor_res$estimate)
  p_values <- c(p_values,cor_res$p.value)
}

names(r_values) <- rownames(bacteria_df)
names(p_values) <- rownames(bacteria_df)
p_adjusted <- p.adjust(p_values,method = "BH")
res <- data.frame(r=r_values,
                  p=p_values,
                  p_adjusted=p_adjusted)

res %>% filter(p_adjusted<0.05)



