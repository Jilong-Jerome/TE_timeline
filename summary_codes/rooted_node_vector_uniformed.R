library(tidyverse)
library(ape)
library(ggpubr)
args <- commandArgs(T)
rooted_tree <- args[1]
node_id <-args[2]
social_sp <- args[3]
subsocial_sp <- args[4]
speciation_thres <- as.numeric(args[5])/100

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

# Replace with the path to your Newick tree file
tree_file_path <- args[1]
node_id <-args[2]
# Read the tree from the file
phylo_tree <- read.tree(file = tree_file_path)
# Assuming phylo_tree is your phylogenetic tree object

# Terminal node IDs (implicitly numbered from 1 to the number of terminal nodes)
node_ids <- 1:length(phylo_tree$tip.label)

# Extract terminal node names
terminal_nodes <- phylo_tree$tip.label

# Species list
species_list <- c("dum", "ten", "sar", "bic", "mim", "lin")

# Function to create species vector
create_species_vector <- function(node_name) {
  parts <- unlist(strsplit(node_name, "\\|")) # Split by '|'
  species_name <- unlist(strsplit(parts[2], "_"))[1] # Split the second part by '_'
  as.integer(species_list %in% species_name) # Create binary vector
}

# Parse species information and create vectors
species_info <- sapply(terminal_nodes, create_species_vector)

# Function to recursively calculate species vector for internal nodes
calculate_species_vector <- function(node, tree) {
  if(node <= length(tree$tip.label)) {
    # If terminal node, return its species vector
    return(species_info[, node])
  } else {
    # Find child nodes
    children <- which(tree$edge[,1] == node)
    child_nodes <- tree$edge[children, 2]
    
    # Recursively calculate for child nodes
    child_vectors <- lapply(child_nodes, function(x) calculate_species_vector(x, tree))
    
    # Sum and average the child vectors
    if(length(child_vectors) == 2) {
      return((child_vectors[[1]] + child_vectors[[2]]) / 2)
    } else {
      # Return NA vector if not exactly two children
      return(rep(NA, length(species_list)))
    }
  }
}

# Calculate species vectors for internal nodes
internal_node_ids <- (length(phylo_tree$tip.label) + 1):NROW(phylo_tree$edge)
internal_node_names <- paste("InternalNode", internal_node_ids)
internal_species_vectors <- sapply(internal_node_ids, function(x) calculate_species_vector(x, phylo_tree), simplify = "array")

# Combine terminal and internal nodes
all_node_ids <- c(node_ids, internal_node_ids)
all_node_names <- c(terminal_nodes, internal_node_names)
all_species_vectors <- cbind(species_info, internal_species_vectors)

# Create data frame
node_data <- data.frame(
  node_id = all_node_ids,
  node_name = all_node_names,
  species_vector = I(split(all_species_vectors, col(all_species_vectors)))
)

# Calculate all pairwise distances between nodes
all_distances <- dist.nodes(phylo_tree)
# Function to calculate the shortest distances and save results
calculate_and_save_shortest_distances <- function(tree, distances) {
    results <- data.frame(Node = character(), Distance = numeric(), stringsAsFactors = FALSE)
    num_tips <- length(tree$tip.label)
    total_nodes <- num_tips + tree$Nnode
    for (internal_node in (num_tips+1):total_nodes) {
        # Ensure the index is within the bounds of the matrix
        if (internal_node <= nrow(distances)) {
            # Distances from this internal node to all terminal nodes
            distances_to_terminals <- distances[internal_node, 1:num_tips]
            shortest_distance <- min(distances_to_terminals)
            results <- rbind(results, data.frame(Node = paste("InternalNode", internal_node), Distance = shortest_distance))
        }
    }
    return(results)
}
# Calculate the shortest distances and save them
distances_df <- calculate_and_save_shortest_distances(phylo_tree, all_distances)
distances_df <- distances_df %>%
      mutate(Distance = ifelse(Distance >0, Distance, 0))
# Find the largest shortest distance
largest_shortest_distance <- max(distances_df$Distance)
print(paste("Largest shortest distance:", largest_shortest_distance))

# Merge the distances_df with node_data
node_data_extended <- merge(node_data, distances_df, by.x = "node_name", by.y = "Node", all.x = TRUE)

# Rename the columns appropriately
names(node_data_extended)[names(node_data_extended) == "Distance"] <- "Shortest_Distance_to_Terminal"

# Replace NA in Shortest_Distance_to_Terminal for terminal nodes with NA_real_
node_data_extended$Shortest_Distance_to_Terminal[is.na(node_data_extended$Shortest_Distance_to_Terminal)] <- NA_real_

# Function to extract terminal node indices for a given species
get_species_indices <- function(tree, species_code) {

# Match species code which comes after a "|" and before an "_"
    species_nodes <- grep(paste0("\\|", species_code, "_"), tree$tip.label)
    return(species_nodes)
}
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

node_data_extended <- node_data_extended%>%
	mutate(uniformed_time = Shortest_Distance_to_Terminal/speciation_scaler)

write_tsv(node_data_extended,paste(node_id,"_",speciation_thres,"_scaled_timeline.tsv",sep=""))


