# Network analysis
setwd("~/IKEM/Projects/PSC_study/")
source("analysis/scripts/custom_functions.R")

# Data import ----
## Czech
path = "data/analysis_ready_data/ikem/"
asv_tab_ikem <- as.data.frame(fread(file.path(path,"asv_table_ikem.csv"),
                                    check.names = FALSE))
taxa_tab_ikem <- as.data.frame(fread(file.path(path,"taxa_table_ikem.csv"),
                                     check.names = FALSE))
metadata_ikem <- as.data.frame(fread(file.path(path,"metadata_ikem.csv"),
                                     check.names = FALSE))

## Norway
path = "data/analysis_ready_data/norway/"
asv_tab_norway <- as.data.frame(fread(file.path(path,"asv_table_norway.csv"),
                                      check.names = FALSE))
taxa_tab_norway <- as.data.frame(fread(file.path(path,"taxa_table_norway.csv"),
                                       check.names = FALSE))
metadata_norway <- as.data.frame(fread(file.path(path,"metadata_norway.csv"),
                                       check.names = FALSE))

## Merging ----
### Terminal ileum
  
ileum_data <- merging_data(asv_tab_1=asv_tab_ikem,
                           asv_tab_2=asv_tab_norway,
                           taxa_tab_1=taxa_tab_ikem,
                           taxa_tab_2=taxa_tab_norway,
                           metadata_1=metadata_ikem,
                           metadata_2=metadata_norway,
                           segment="TI",Q="Q1")

ileum_asv_tab <- ileum_data[[1]]
ileum_taxa_tab <- ileum_data[[2]]
ileum_metadata <- ileum_data[[3]]

### Colon
  
colon_data <- merging_data(asv_tab_1=asv_tab_ikem,
                           asv_tab_2=asv_tab_norway,
                           taxa_tab_1=taxa_tab_ikem,
                           taxa_tab_2=taxa_tab_norway,
                           metadata_1=metadata_ikem,
                           metadata_2=metadata_norway,
                           segment="colon",Q="Q1")

colon_asv_tab <- colon_data[[1]]
colon_taxa_tab <- colon_data[[2]]
colon_metadata <- colon_data[[3]]

# Aggregation, filtering
level="genus"
## Ileum
genus_data <- aggregate_taxa(ileum_asv_tab,
                             ileum_taxa_tab,
                             taxonomic_level=level,
                             names=TRUE)
# Filtration
filt_data <- filtering_steps(genus_data[[1]],
                             genus_data[[2]],
                             ileum_metadata,
                             seq_depth_threshold=10000)

filt_ileum_genus_tab <- filt_data[[1]]
filt_ileum_genus_taxa <- filt_data[[2]]
filt_ileum_metadata <- filt_data[[3]]

## Colon
genus_data <- aggregate_taxa(colon_asv_tab,
                             colon_taxa_tab,
                             taxonomic_level=level,
                             names=TRUE)

filt_data <- filtering_steps(genus_data[[1]],
                             genus_data[[2]],
                             colon_metadata,
                             seq_depth_threshold=10000)

filt_colon_genus_tab <- filt_data[[1]]
filt_colon_genus_taxa <- filt_data[[2]]
filt_colon_genus_metadata <- filt_data[[3]]

# Network analysis ----

## Terminal ileum ----
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
psc_taxa <- read.xlsx("analysis/results/Q1/univariate_analysis/psc_effect_terminal_ileum.xlsx")
psc_taxa <- psc_taxa$SeqID

# Check if rownames and colnames are in names_vector
rownames_in_vector <- rownames(similarity_matrix) %in% psc_taxa
colnames_in_vector <- colnames(similarity_matrix) %in% psc_taxa

# Create the boolean matrix
boolean_matrix <- outer(rownames_in_vector, colnames_in_vector, `&`)

# Assign rownames and colnames back to the boolean matrix
dimnames(boolean_matrix) <- dimnames(similarity_matrix)

# View the boolean matrix
View(boolean_matrix)

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

# Add attributes for visualization
g <- graph_from_data_frame(edge_list, directed = FALSE)
V(g)$size <- 15
V(g)$color <- "skyblue"
E(g)$width <- abs(edge_list$Weight) * 5  # Scale edge width by weight
E(g)$color <- ifelse(edge_list$Weight > 0, "#3f00fb", "#f32b1d")  # Green for positive, red for negative

# Save edge list
write_csv(edge_list, "analysis/results/Q1/network/terminal_ileum_taxon_network_edges.csv")

# Save node attributes
psc_effect_ileum_genus <- read.xlsx("analysis/results/Q1/univariate_analysis/psc_effect_terminal_ileum.xlsx")
increased <- psc_effect_ileum_genus$SeqID[psc_effect_ileum_genus$log2FoldChange >0]
decreased <- psc_effect_ileum_genus$SeqID[psc_effect_ileum_genus$log2FoldChange <0]
nodes <- V(g)$name
sizes <- rep(10,length(V(g)))
colors <- rep("neutral",length(V(g)))
colors[nodes %in% increased] <- "increased"
colors[nodes %in% decreased] <- "decreased"
              
node_attributes <- data.frame(
  Node = nodes,
  Size = sizes, 
  Color = colors # increased or deacreasing 
)

write_csv(node_attributes, "analysis/results/Q1/network/terminal_ileum_taxon_network_nodes.csv")


## Colon ----
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
psc_taxa <- read.xlsx("analysis/results/Q1/univariate_analysis/psc_effect_colon.xlsx")
psc_taxa <- psc_taxa$SeqID

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

# Add attributes for visualization
V(g)$size <- 15
V(g)$color <- "skyblue"
E(g)$width <- abs(edge_list$Weight) * 5  # Scale edge width by weight
E(g)$color <- ifelse(edge_list$Weight > 0, "#3f00fb", "#f32b1d")  # Green for positive, red for negative

# Save edge list
write_csv(edge_list, "analysis/results/Q1/network/colon_taxon_network_edges.csv")

# Save node attributes
psc_effect_colon_genus <- read.xlsx("analysis/results/Q1/univariate_analysis/psc_effect_colon.xlsx")
increased <- psc_effect_colon_genus$SeqID[psc_effect_colon_genus$log2FoldChange >0]
decreased <- psc_effect_colon_genus$SeqID[psc_effect_colon_genus$log2FoldChange <0]
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

write_csv(node_attributes, "analysis/results/Q1/network/colon_taxon_network_nodes.csv")
