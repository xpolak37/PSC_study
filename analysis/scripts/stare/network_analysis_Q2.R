# network analysis
setwd("~/IKEM/norsko/data analysis")

# Q2 genus level

library(ccrepe)
library(igraph)
library(readr)

# Example data matrix (samples x taxa)
data_matrix <- filt_ileum_genus_tab %>% column_to_rownames("SeqID") %>% t() %>% as.data.frame()
data_matrix <- t(apply(data_matrix,1, function(x) x/sum(x)))

# Run CCREPE
ccrepe_results_ileum <- ccrepe(data_matrix)

# Extract similarity and p-value matrices
similarity_matrix <- ccrepe_results_ileum$sim.score  # Association strengths
p_value_matrix <- ccrepe_results_ileum$p.values     # P-values

similarity_matrix[is.na(similarity_matrix)] <- 0
p_value_matrix[is.na(p_value_matrix)] <- 1

# Define thresholds for filtering
similarity_threshold <- 0.3  # Minimum association strength
p_value_threshold <- 1.02   # Maximum p-value for significance
psc_taxa <- read.xlsx("results/Q2/rPSC_effect_terminal_ileum.xlsx",sheet = "Ileum Genus")
psc_taxa <- psc_taxa$ASV


# Check if rownames and colnames are in names_vector
rownames_in_vector <- rownames(similarity_matrix) %in% psc_taxa
colnames_in_vector <- colnames(similarity_matrix) %in% psc_taxa

# Create the boolean matrix
boolean_matrix <- outer(rownames_in_vector, colnames_in_vector, `&`)

# Assign rownames and colnames back to the boolean matrix
dimnames(boolean_matrix) <- dimnames(similarity_matrix)

# View the boolean matrix
boolean_matrix


# Filter associations
significant_edges <- which(
  p_value_matrix < p_value_threshold, 
  arr.ind = TRUE
)

# Create edge list for significant associations
edge_list <- data.frame(
  Source = colnames(similarity_matrix)[significant_edges[, 1]],
  Target = rownames(similarity_matrix)[significant_edges[, 2]],
  Weight = similarity_matrix[significant_edges],
  width = abs(similarity_matrix[significant_edges]),
  color = similarity_matrix[significant_edges]>0,
  significant=p_value_matrix[significant_edges]<0.05
)

# Remove self-loops (if any)
edge_list <- edge_list[edge_list$Source != edge_list$Target, ]

# Add canonical "key" columns to deduplicate
edge_list_canonical <- edge_list
edge_list_canonical$key1 <- pmin(edge_list$Source, edge_list$Target)
edge_list_canonical$key2 <- pmax(edge_list$Source, edge_list$Target)

# Remove duplicate rows based on keys
edge_list_nonredundant <- edge_list_canonical[!duplicated(edge_list_canonical[c("key1", "key2")]), ]

# Drop the temporary key columns
edge_list_nonredundant <- edge_list_nonredundant[, !(names(edge_list_nonredundant) %in% c("key1", "key2"))]

edge_list <- edge_list_nonredundant

# Add attributes for visualization (optional)
g <- graph_from_data_frame(edge_list, directed = FALSE)
V(g)$size <- 15
V(g)$color <- "skyblue"
E(g)$width <- abs(edge_list$Weight) * 5  # Scale edge width by weight
E(g)$color <- ifelse(edge_list$Weight > 0, "#3f00fb", "#f32b1d")  # Green for positive, red for negative

# Save edge list
write_csv(edge_list, "results/Q2/network/Genus terminal ileum/taxon_network_edges.csv")

# Save node attributes
psc_effect_ileum_genus <- read.xlsx("results/Q2/rPSC_effect_terminal_ileum.xlsx",sheet = "Ileum Genus")
increased <- psc_effect_ileum_genus$ASV[psc_effect_ileum_genus$log2FoldChange >0]
decreased <- psc_effect_ileum_genus$ASV[psc_effect_ileum_genus$log2FoldChange <0]
nodes <- V(g)$name
sizes <- rep(10,length(V(g)))
colors <- rep("neutral",length(V(g)))
colors[nodes %in% increased] <- "increased"
colors[nodes %in% decreased] <- "decreased"

node_attributes <- data.frame(
  Node = nodes,
  Size = sizes, 
  Color = colors # podla increased or deacreasing LOGFOLDCHANGE ZNAMIENKO
)
write_csv(node_attributes, "results/Q2/network/Genus terminal ileum/taxon_network_nodes.csv")


# COLON
# Example data matrix (samples x taxa)
data_matrix <- filt_colon_genus_tab %>% column_to_rownames("SeqID") %>% t() %>% as.data.frame()
data_matrix <- t(apply(data_matrix,1, function(x) x/sum(x)))

# Run CCREPE
ccrepe_results_colon <- ccrepe(data_matrix)

# Extract similarity and p-value matrices
similarity_matrix <- ccrepe_results_colon$sim.score  # Association strengths
p_value_matrix <- ccrepe_results_colon$p.values     # P-values


# Define thresholds for filtering
similarity_threshold <- 0.3  # Minimum association strength
p_value_threshold <- 1.02  # Maximum p-value for significance
psc_taxa <- read.xlsx("results/Q2/rPSC_effect_colon.xlsx",sheet = "colon Genus")
psc_taxa <- psc_taxa$ASV


# Filter associations
significant_edges <- which(
  p_value_matrix < p_value_threshold, 
  arr.ind = TRUE
)

# Create edge list for significant associations
edge_list <- data.frame(
  Source = colnames(similarity_matrix)[significant_edges[, 1]],
  Target = rownames(similarity_matrix)[significant_edges[, 2]],
  Weight = similarity_matrix[significant_edges],
  width = abs(similarity_matrix[significant_edges]),
  color = similarity_matrix[significant_edges]>0,
  significant=p_value_matrix[significant_edges]<0.05
)


# Remove self-loops (if any)
edge_list <- edge_list[edge_list$Source != edge_list$Target, ]
# Add canonical "key" columns to deduplicate
edge_list_canonical <- edge_list
edge_list_canonical$key1 <- pmin(edge_list$Source, edge_list$Target)
edge_list_canonical$key2 <- pmax(edge_list$Source, edge_list$Target)

# Remove duplicate rows based on keys
edge_list_nonredundant <- edge_list_canonical[!duplicated(edge_list_canonical[c("key1", "key2")]), ]

# Drop the temporary key columns
edge_list_nonredundant <- edge_list_nonredundant[, !(names(edge_list_nonredundant) %in% c("key1", "key2"))]

edge_list <- edge_list_nonredundant

g <- graph_from_data_frame(edge_list, directed = FALSE)
# Add attributes for visualization (optional)
V(g)$size <- 15
V(g)$color <- "skyblue"
E(g)$width <- abs(edge_list$Weight) * 5  # Scale edge width by weight
E(g)$color <- ifelse(edge_list$Weight > 0, "#3f00fb", "#f32b1d")  # Green for positive, red for negative

# Save edge list
write_csv(edge_list, "results/Q2/network/Genus colon/taxon_network_edges.csv")

# Save node attributes
psc_effect_colon_genus <- read.xlsx("results/Q2/rPSC_effect_colon.xlsx",sheet = "colon Genus")
increased <- psc_effect_colon_genus$ASV[psc_effect_colon_genus$log2FoldChange >0]
decreased <- psc_effect_colon_genus$ASV[psc_effect_colon_genus$log2FoldChange <0]
nodes <- V(g)$name
sizes <- rep(10,length(V(g)))
colors <- rep("neutral",length(V(g)))
colors[nodes %in% increased] <- "increased"
colors[nodes %in% decreased] <- "decreased"

node_attributes <- data.frame(
  Node = nodes,
  Size = sizes, 
  Color = colors # podla increased or deacreasing LOGFOLDCHANGE ZNAMIENKO
)
write_csv(node_attributes, "results/Q2/network/Genus colon/taxon_network_nodes.csv")
