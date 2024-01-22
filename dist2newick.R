library(ape)
parse_root_set <- function(file_path) {
  # Read the file
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Check if the file has the expected number of rows and columns
  if (nrow(data) != 1 || ncol(data) != 5) {
    stop("The file format is not as expected.")
  }

  # Extracting individual parameters
  node_id <- data$Node.ID
  focus_species <- data$Focus.Species
  closest_species <- data$Closest.Species
  selected_sequence <- data$Selected.Sequence
  root_situation <- data$Root.Situation

  # Return as a list
  return(list(node_id = node_id, focus_species = focus_species, closest_species = closest_species, selected_sequence = selected_sequence, root_situation = root_situation))
}


args <- commandArgs(T)
# Read the distance matrix
# Replace 'your_file.csv' with the path to your distance matrix file
# Adjust the read function and its parameters as per your file format
dist_matrix <- as.matrix(read.csv(args[1], row.names = 1))
root_res <- parse_root_set(args[2])
node <- root_res$node_id
sp <- root_res$focus_species
out_sp <- root_res$closest_species
out_seq <- root_res$selected_sequence
root_info <- root_res$root_situation

# Build the phylogenetic tree using Neighbor-Joining method
tree <- nj(dist_matrix)
# Check if the tree is rooted
is_rooted <- !is.null(tree$root.edge) && length(tree$root.edge) > 0

# Print the result
if (is_rooted) {
  cat("The tree is rooted.\n")
} else {
  cat("The tree is unrooted.\n")
}

# Plot and save the unrooted tree
pdf(paste(sp,"_",gsub("_R_","", node),"_unrooted.pdf",sep=""))  # Save the plot to a PDF file
plot(tree,show.tip.label = FALSE,type = "unrooted")
dev.off()  # Close the PDF device
write.tree(tree, file = paste(sp,"_",gsub("_R_","", node),"_unrooted.newick",sep=""))

# Specify the outgroup label
# Extract tip labels
tip_labels <- tree$tip.label
# Find labels containing "mim"
labels_hit <- tip_labels[grep(out_seq, tip_labels)]
outgroup_label <- labels_hit[1]

# Check if the outgroup label exists in the tree
if (!outgroup_label %in% tree$tip.label) {
    stop("Specified outgroup label not found in the tree tip labels.")
}

# Root the tree with the specified outgroup
rooted_tree <- root(tree, outgroup = outgroup_label, resolve.root = TRUE)
# Plot and save the unrooted tree
pdf(paste(sp,"_",gsub("_R_","", node),"_rooted.pdf",sep=""))  # Save the plot to a PDF file
plot(rooted_tree,show.tip.label = FALSE)
dev.off()  # Close the PDF device
write.tree(rooted_tree, file = paste(sp,"_",gsub("_R_","", node),"_rooted.newick",sep=""))
