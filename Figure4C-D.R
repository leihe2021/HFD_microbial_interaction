#----------------------------------
# figure 4C alpha diversity ----
#----------------------------------
rm(list = ls())
library(tidyverse)
library(microeco)
library(magrittr)
load("otu_grp_taxa.Rda")
dataset <- microtable$new(sample_table = group,
                          otu_table = otu, 
                          tax_table = taxa)
dataset$sample_table$group <- factor(dataset$sample_table$group, 
                                     levels = c("Control", "HFD", "FO"))
t1 <- trans_alpha$new(dataset = dataset, group = "group")
input <- t1$data_alpha %>% 
  filter(Measure == "Chao1") %>% 
  select(group, sample_id, Measure, Value)

ggplot(input, aes(x = group, y = Value, colour = group)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(size = 3, width = .2) +
  scale_color_manual(values = c("skyblue","red", "purple")) +
  theme_classic() +
  theme(legend.position = "none")
#----------------------------------
# figure 4D PCo plot
# calculate beta-diversity
#----------------------------------
dataset$cal_betadiv()
t1 <- trans_beta$new(dataset = dataset, group = "group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
t1$plot_ordination(plot_color = "group", plot_type = c("point", "ellipse")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("skyblue","red", "purple"))
