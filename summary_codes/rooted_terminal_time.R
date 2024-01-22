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

# Replace with the path to your Newick tree file
tree_file_path <- "N1168_blast_rooted.newick"
node_id <- "N1168"
# Read the tree from the file
phylo_tree <- read.tree(file = tree_file_path)
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

# View the resulting dataframe
print(max(df$TerminalBranchLength))
print(min(df$TerminalBranchLength))

# Determine the range (0 to maximum value of ee) and create 100 intervals
max_length <- max(branch)
breaks <- seq(0, max_length, length.out = 101)
intervals <- cut(branch, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# Calculate interval start and end points
interval_starts <- head(breaks, -1)
interval_ends <- tail(breaks, -1)

# Add interval information to the dataframe
df$IntervalID <- intervals
df$IntervalStart <- interval_starts[intervals]
df$IntervalEnd <- interval_ends[intervals]

df <- df %>%
  mutate(Chromosome = sapply(strsplit(NodeName, "\\|"), function(x) x[2]),
         Species = sapply(strsplit(Chromosome, "_"), function(x) x[1]),
         PositionStart = as.numeric(sapply(strsplit(NodeName, "\\|"), function(x) x[3])),
         PositionEnd = as.numeric(sapply(strsplit(NodeName, "\\|"), function(x) x[4])),
         MeanPosition = (PositionStart + PositionEnd) / 2)

# Remove temporary columns if you don't need them
df <- df %>%
  select(-PositionStart, -PositionEnd)
# View the resulting dataframe
print(df)

# Calculate all pairwise distances between nodes
all_distances <- dist.nodes(phylo_tree)
# Function to extract terminal node indices for a given species
get_species_indices <- function(tree, species_code) {
    # Match species code which comes after a "|" and before an "_"
    species_nodes <- grep(paste0("\\|", species_code, "_"), tree$tip.label)
    return(species_nodes)
}
# Extract indices of terminal nodes for species "dum" and "ten"
dum_indices <- get_species_indices(phylo_tree, "dum")
ten_indices <- get_species_indices(phylo_tree, "ten")
# Extract pairwise distances between "dum" and "ten" terminal nodes
# This will only include inter-species distances as we are selecting between two different species
pairwise_distances <- all_distances[dum_indices, ten_indices]
# Convert the matrix to a vector, excluding diagonal elements (intra-species distances)
pairwise_distances_vector <- pairwise_distances[upper.tri(pairwise_distances)]

# Adjust pairwise distances by dividing by 2
adjusted_pairwise_distances_vector <- pairwise_distances_vector / 2
sp_low <- quantile(adjusted_pairwise_distances_vector, 0.01)
sp_high <- quantile(adjusted_pairwise_distances_vector, 0.05)

# Plot single figure
fai_data <- read.table("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/all_6_genome.fa.fai",header = FALSE, col.names = c("Chromosome", "Length","X1","X2","X3"))
# Extract species code from chromosome ID and add it to the fai_data
fai_data$Species <- gsub("_.*", "", fai_data$Chromosome)

temp <- df %>%
	group_by(Species)%>%
	summarize(n = n())
max_count <- max(temp$n)
count_dt_init = data.frame(
    Species = c("dum","ten","sar","bic","mim","lin"),
    total = c(0,0,0,0,0,0)
)
for (i in 1:100){

    df_toplot <- df %>% filter(IntervalID >= i)
    current_time <- min(df_toplot$IntervalStart)
    oldest_time <- max_length
    timeplot <-
	    ggplot()+
	    geom_segment(aes(x = 0, xend = oldest_time,y=1,yend=1),size = 2,alpha = 0, color = "grey" )+
	    geom_segment(aes(x = current_time, xend = oldest_time,y=1,yend=1),size = 2,alpha = 0.8, color = "#657CD3" )+
	    geom_segment(aes(x = sp_low,xend=sp_high,y=1,yend=1),size = 2, alpha = 1, color = "red")+
	    geom_text(aes(x = (sp_low+sp_high)/2, y = 1.5, label = "dum ten divergence"))+
	    geom_segment(aes(x = (sp_low+sp_high)/2, xend = (sp_low+sp_high)/2, y = 1.4,yend = 1.1),arrow = arrow(type = "open"))+
	    xlab("Relative Time")+
	    ylab("")+
	    theme(axis.ticks.y = element_blank(),
               axis.text.y = element_blank())

    p <- ggplot()+
         geom_segment(data=fai_data,mapping = aes(x=1, xend = Length, y = Chromosome, yend = Chromosome),color = "grey", alpha = 0.6) +
         geom_point(data = df_toplot%>% filter(!is.na(MeanPosition)), mapping = aes(x = MeanPosition, y = Chromosome,color = Species))+
         facet_wrap(~ factor(Species,levels = c("dum","sar","mim","ten","bic","lin")),scales = "free")+
         scale_colour_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey"))+
        scale_alpha_continuous(range = c(0, 1))+
        ggtitle(node_id) +
        xlab("Position") +
        ylab("Chromosome ID")+
	theme(legend.position = "none")
    sum_df <- df_toplot%>% filter(!is.na(MeanPosition))%>%
	    group_by(Species)%>%
	    summarize(total = n())
    sum_count <- ggplot()+
	    geom_col(data = count_dt_init, mapping = aes(x= factor(Species, levels = c("dum","ten","sar","bic","mim","lin")),y = total ,color = Species,fill = Species))+
	    geom_col(data = sum_df, mapping = aes(x= factor(Species,levels = c("dum","ten","sar","bic","mim","lin")),y = total ,color = Species,fill = Species))+
	    ylim(0,max_count)+
	    xlab("Species")+
	    scale_colour_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey"))+
	    scale_fill_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey"))
    combo_1 <- ggarrange(timeplot,p,ncol=1, nrow= 2, heights = c(0.2,1))
    combo_2 <- ggarrange(combo_1,sum_count,ncol=2, nrow= 1, heights = c(1,0.2))
    # Save the plot as a PNG file
    png_filename <- paste0("Interval_",i, ".png")
    ggsave(png_filename, plot = combo_2, width = 10, height = 10)
}
