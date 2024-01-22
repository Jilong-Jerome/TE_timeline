library(tidyverse)
library(ggpubr)

args <- commandArgs(T)
node <- args[1]
sp <- args[2]
sister_sp <- args[3]
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

count_df <- read_tsv(paste(node,"_interval_sum.tsv",sep=""))%>%
  pivot_longer(cols = c("dum", "ten", "sar", "bic", "mim", "lin"),
               names_to = "Species",
               values_to = "Count") %>%
  mutate(Species = factor(Species))%>%
  filter(Count > 0)

count_exp <- ggplot(count_df, aes(x = Interval_ID, y = Count, color = Species)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 26, linetype = "dashed", color = "red")+
  scale_colour_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey")) +
  ylab("Expansion_Count") +
  xlab("")+
  ggtitle(node)

fraction_df <- read_tsv(paste(node,"_interval_sum.tsv",sep=""))%>%
  mutate(total = dum + ten + sar+ bic+ mim + lin)%>%
  mutate(dum = dum/total,
         ten = ten/total,
         sar = sar/total,
         bic = bic/total,
         mim = mim/total,
         lin = lin/total)%>%
  pivot_longer(cols = c("dum", "ten", "sar", "bic", "mim", "lin"),
               names_to = "Species",
               values_to = "Fraction") %>%
  mutate(Species = factor(Species,levels =c("dum", "ten", "sar", "bic", "mim", "lin")))

fraction_exp <- fraction_df%>%
  ggplot(aes(x = Interval_ID, y = Fraction, color = Species,fill = Species)) +
  geom_col(alpha = 0.7) +
  scale_colour_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey")) +
  scale_fill_manual(name = "Species",
                      values = c("dum" = "#F94C10",
                                 "ten" = "#362FD9",
                                 "sar" = "#C70039",
                                 "bic" = "#2E97A7",
                                 "mim" = "#900C3F",
                                 "lin" = "grey")) +
  ylab("Expansion_fraction") +
  xlab("Time_interval (Left:Current, Right:Past)")
  #ggtitle("N1168:Retrovirus-like")

combo_fig <- ggarrange(count_exp,fraction_exp,ncol =  1)
ggsave(paste(sp,"_",node,"_timeline.png",sep=""),combo_fig)

pair_d <- read_tsv(paste(sp,"_",sister_sp,"_pair_dist.tsv",sep=""))

div_fig <- pair_d %>%
  ggplot()+
  geom_segment(aes(x=0,xend = 1,y = 0, yend = pair_d),alpha = 0.005,color = "#3258a8")+
  geom_segment(aes(x=1,xend = 2,y = pair_d, yend = 0),alpha = 0.005, color = "#3258a8")+
  geom_hline(yintercept = quantile(pair_d$pair_d,0.01),color ="red")+
  ylab("distance")+
  ggtitle(paste(node,sp,"-",sister_sp,"pairwise distance"))+
  theme(
    axis.title.x = element_blank(),  # Hide x-axis label
    axis.text.x = element_blank(),   # Hide x-axis text/ticks
    axis.ticks.x = element_blank()   # Hide x-axis ticks
)
ggsave(paste(node,"_sp_div.png",sep=""),div_fig)


