library(tidyverse)
library(network)
library(intergraph)
library(GGally)
library(colorRamps)
library(igraph)
library(intergraph)
library(ggpubr)

# read structural varaints dataset
read_sp_sv<- function(species_sv) {
  sp_sv_df <- read_tsv(species_sv,col_names = c("sp","chrom","start","end","len","type","ID","seqID"))
  return(sp_sv_df)
}

# draw SV length distribtuions
sv_pair_plot <- function(social_sv,subsocial_sv,social,subsocial) {
  social_df <- read_sp_sv(social_sv)%>%filter(len >=100)%>%filter(len <=10000)
  subsocial_df <- read_sp_sv(subsocial_sv)%>%filter(len >=100)%>%filter(len <=10000)
  p <- ggplot() + 
    geom_histogram(data = social_df, aes_string(x = "len", fill = shQuote(social)), alpha = 0.4, bins = 200) +
    geom_histogram(data = subsocial_df, aes_string(x = "len", fill = shQuote(subsocial)), alpha = 0.4, bins = 200) +
    scale_fill_manual(values = setNames(c("red", "blue"), c(social,subsocial)),
                      name = "Species",
                      breaks = c(social,subsocial))+
    xlab("length")+
    labs(title = "Overlaying Histograms from Structural Variants Length by type")+
    facet_wrap(~type)
  return(p)
}

# draw SV graph
sv_graph <- function(rds,social,subsocial){
# color panel
color_ramp <- colorRampPalette(c("#0000FF", "grey", "#FF0000"))
colors.bitrait = color_ramp(11)
# load graph data
nets <- readRDS(rds)
g.part <- nets$net
node_df <- nets$node_df
seq_df <- nets$seq_df
sp_f <- nets$species_f
temp_g <- asIgraph(g.part)
# Get vertex ids with node_size > 10
v_ids <- V(temp_g)[V(temp_g)$node_size >= 10]
# Get Subgraph
sub_g <- induced_subgraph(temp_g, v_ids)
sub_net <- asNetwork(sub_g)
sub.graph.names = network.vertex.names(sub_net)
set.seed(239)
p.refined <- ggnet2(sub_net, size= "node_size",
       color=  ifelse(paste(">",as.character(as.integer(sp_f[sub.graph.names]/0.1)*10),"%",sep="")==">100%",">90%",paste(">",as.character(as.integer(V(sub_g)$f_dum/0.1)*10),"%",sep="")),
       arrow.gap = 0.01, arrow.size = 2,
       palette = c(">0%" = colors.bitrait[1],
                        ">10%" = colors.bitrait[2],
                        ">20%" = colors.bitrait[3],
                        ">30%" = colors.bitrait[4],
                        ">40%" = colors.bitrait[5],
                        ">50%" = colors.bitrait[7],
                        ">60%" = colors.bitrait[8],
                        ">70%" = colors.bitrait[9],
                        ">80%" = colors.bitrait[10],
                        ">90%" = colors.bitrait[11]),
       label = F)+
   guides(size = F)+
  labs(title = paste("SV graph of ",social," & ",subsocial, "(node size >= 10)",sep=""))
return(p.refined)
}

# count number of sequences and nodes for each species in a graph
svg_number <- function(rds){
nets <- readRDS(rds)
g.part <- nets$net
node_df <- nets$node_df
seq_df <- nets$seq_df
sp_f <- nets$species_f
large_node<-
  node_df%>%
  filter(n >= 10)
svg_num <- seq_df%>%
  filter(node %in% large_node$node)%>%
  separate(name,into = c("sp","chr","pos","type","status","id","size"),remove = F,sep = "\\|")%>%
  group_by(sp)%>%
  summarise( n_seq= n(),
             n_node = length(unique(node)))%>%
  mutate(total_node = NROW(large_node))
return(svg_num)
}

# historgram of node number based on fraction of social species sequence
svg_hist <- function(rds,social,subsocial){
node_summary <- node_sum(rds,social,subsocial)
p_hist <- node_summary%>%
  ggplot(aes(x=dum_f))+
  geom_histogram(bins = 100)+
  xlab(paste("Fraction of social sequences"))+
  labs(title = paste("Histogram of social fraction in nodes for ",social, " & ",subsocial,sep=""))+
  ylab("number of nodes")
return(p_hist)
}

# Unique node size jitter plot
unique_jitter<-function(rds,social,subsocial){
node_summary <- node_sum(rds,social,subsocial)
p_jitter<- node_summary%>%
    mutate(species = ifelse(dum_f >= 1,social,ifelse(dum_f <= (1-1), subsocial,"mix")))%>%
  filter(species %in% c(social,subsocial))%>%
  group_by(species)%>%
  mutate(species = factor(species, levels = c(social,subsocial)))%>%
  ggplot()+
  geom_jitter(aes(x = species,y=size,color = species))+
  ylab("node size")+
  scale_color_manual(values = setNames(c("red", "blue"), c(social,subsocial)),
                      name = "Species",
                      breaks = c(social,subsocial))
return(p_jitter)
}

#node_summary
node_sum<- function(rds,social,subsocial) {
  nets <- readRDS(rds)
g.part <- nets$net
node_df <- nets$node_df
seq_df <- nets$seq_df
sp_f <- nets$species_f
large_node<-
  node_df%>%
  filter(n >= 10)
node_summary <- seq_df%>%
  filter(node %in% large_node$node)%>%
  separate(name,into = c("sp","chr","pos","type","status","id","size"),remove = F,sep = "\\|")%>%
  #filter(type %in% c("DUP","INVDP"))%>%
  group_by(node)%>%
  summarise(
    count_dum = sum(sp == social),
    count_ten = sum(sp == subsocial),
    chrom_range = n_distinct(chr),
    max_len = max(as.numeric(size))
  )%>%
  mutate(size = (count_dum + count_ten),
         dum_f = count_dum / size,
         dum_f_se = sqrt(dum_f*(1-dum_f)/size),
         multi_chrom = ifelse(chrom_range > 1, TRUE,FALSE))
return(node_summary)
}
# sequence divergence in a list of nodes
# sequence divergence in a list of nodes
node_div<-function(div_df,node_list,social,subsocial){
div_p <- div_df%>%
  filter(seq1_name != seq2_name)%>%
  filter(node %in% node_list)%>%
  group_by(`species-pair`)%>%
  #filter(`species-pair` %in% c(paste(social,social,sep = "_"),paste(subsocial,subsocial,sep = "_")))%>%
  ggplot()+
  geom_histogram(aes(x=average_identity,fill= `species-pair`),bins = 200)+
  facet_wrap(~`species-pair`,nrow = 3)+
  scale_fill_manual(values = setNames(c("red", "blue"), c(paste(social,social,sep = "_"),paste(subsocial,subsocial,sep = "_"))),
                      name = "Species",
                      breaks = c(paste(social,social,sep = "_"),paste(subsocial,subsocial,sep = "_")))

# if there is cross species comparision, spciation is defined as 95% quantile
species_pair <- sort(c(social,subsocial))
cross_name <- paste(species_pair[1],species_pair[2],sep="_")
temp_df <- div_df%>%
  filter(node %in% node_list)
if (cross_name %in% unique(temp_df$`species-pair`)){
  div_vec <- temp_df%>%filter(`species-pair` == cross_name)
  t_split <- quantile(div_vec$average_identity,0.95)
} else {
  t_split <- 0
}
  
total_n_df <- div_df%>%
  filter(seq1_name != seq2_name)%>%
  filter(node %in% node_list)%>%
  group_by(`species-pair`)%>%
  summarise(
    n_total = n()
  )%>%
  mutate(pair = `species-pair`)%>%
  select(pair,n_total)
post_n_df <- div_df%>%
  filter(seq1_name != seq2_name)%>%
  filter(node %in% node_list)%>%
  filter(average_identity > t_split)%>%
  group_by(`species-pair`)%>%
  summarise(
    n_post = n()
  )%>%
  mutate(pair = `species-pair`)%>%
  select(pair,n_post)

merge_df <- left_join(post_n_df,total_n_df,by = "pair")
if (paste(social,social,sep="_") %in% merge_df$pair){
  social_post_n <- merge_df%>%filter(pair == paste(social,social,sep="_")) %>%select(n_post)%>%pull()
  social_total_n <- merge_df%>%filter(pair == paste(social,social,sep="_")) %>%select(n_total)%>%pull()
}else{
  social_post_n <- 0
  social_total_n <- 0
}
if (paste(subsocial,subsocial,sep="_") %in% merge_df$pair){
  subsocial_post_n <- merge_df%>%filter(pair == paste(subsocial,subsocial,sep="_")) %>%select(n_post)%>%pull()
  subsocial_total_n <- merge_df%>%filter(pair == paste(subsocial,subsocial,sep="_")) %>%select(n_total)%>%pull()
}else{
  subsocial_post_n <- 0
  subsocial_total_n <- 0
}
df_out <- data.frame(
  node = node_list[1],
  t_split = t_split,
  social_post_n = social_post_n,
  social_total_n = social_total_n,
  social_post_f = social_post_n/social_total_n,
  subsocial_post_n = subsocial_post_n,
  subsocial_total_n = subsocial_total_n,
  subsocial_post_f = subsocial_post_n/subsocial_total_n
)
return(list(p=div_p,ts = t_split,df = df_out))
}

# node seq names
node_seq_name<- function(rds,social,subsocial){
  nets <- readRDS(rds)
  seq_df <- nets$seq_df
  node_sum_df <- node_sum(rds,social,subsocial)
  uniq_node_df <- node_sum_df%>%
  filter(dum_f == 1|dum_f == 0)%>%
  filter(size >= 10)
  seq_out_df <- seq_df%>%
  filter(node %in% uniq_node_df$node)
  return(seq_out_df)
}
