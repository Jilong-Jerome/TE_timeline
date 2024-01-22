library(tidyverse)
library(ape)
library(ggpubr)
# Create a custom theme
my_theme <- function() {
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    legend.position = "right"
  )
}
theme_set(my_theme())
args <- commandArgs(T)
rooted_tree <- args[1]
node_id <-args[2]
social_sp <- args[3]
subsocial_sp <- args[4]
speciation_thres <- as.numeric(args[5])/100

# Read the tree from the file
phylo_tree <- read.tree(file = rooted_tree)
# Identify terminal branches
branchSel <- phylo_tree$edge[,2] < (Ntip(phylo_tree) + 1)

# Extract lengths of terminal branches
branch = phylo_tree$edge.length[branchSel]

# Extract names of terminal nodes
node_names <- phylo_tree$tip.label
print(node_names)
# Creating the dataframe
df <- data.frame(NodeName = node_names, TerminalBranchLength = branch)%>%
	mutate(TerminalBranchLength = ifelse(TerminalBranchLength > 0,TerminalBranchLength,0))

# Function to extract terminal node indices for a given species
get_species_indices <- function(tree, species_code) {

# Match species code which comes after a "|" and before an "_"
    species_nodes <- grep(paste0("\\|", species_code, "_"), tree$tip.label)
    return(species_nodes)
}

# Calculate all pairwise distances between nodes
all_distances <- dist.nodes(phylo_tree)
# Extract indices of terminal nodes for species "dum" and "ten"
dum_indices <- get_species_indices(phylo_tree, social_sp)
ten_indices <- get_species_indices(phylo_tree, subsocial_sp)
# Extract pairwise distances between "dum" and "ten" terminal nodes
# This will only include inter-species distances as we are selecting between two different species
pairwise_distances <- all_distances[dum_indices, ten_indices]

# Convert the matrix to a vector, excluding diagonal elements (intra-species distances)
pairwise_distances_vector <- pairwise_distances[upper.tri(pairwise_distances)]

# Adjust pairwise distances by dividing by 2
adjusted_pairwise_distances_vector <- pairwise_distances_vector / 2
write_tsv(data.frame(pair_d = adjusted_pairwise_distances_vector),paste(social_sp,"_",subsocial_sp,"_pair_dist.tsv",sep=""))
# Calculate 1% and 5% quantiles
speciation_scaler <- quantile(adjusted_pairwise_distances_vector, probs = speciation_thres)

node_data_extended <- df %>%
	mutate(node = node_id)%>%
	mutate(speciation_scaler = speciation_scaler)%>%
        mutate(uniformed_time = TerminalBranchLength/speciation_scaler)

write_tsv(node_data_extended,paste(node_id,"_",speciation_thres,"_scaled_timeline_terminal_branch.tsv",sep=""))
