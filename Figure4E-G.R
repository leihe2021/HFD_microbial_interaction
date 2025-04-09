#------------------------------------------------------------------
# figure 4E barplot for Relative Abundance of bacterial phylum----
#------------------------------------------------------------------
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

t1 <- trans_abund$new(dataset = dataset, taxrank = "phylum", ntaxa = 5)
dim(t1$data_abund)
cols <- c("#FF8381","#923094","#0093C6","#B46C60","#00D472")
t1$plot_bar(color_values = cols, others_color = "grey80", facet = "group", 
            xtext_keep = F, legend_text_italic = F)
#------------------------------------------------------------------
# figure 4F and G barplot for Relative Abundance of bacterial phylum----
#------------------------------------------------------------------
t1 <- trans_diff$new(dataset = dataset, 
                     method = "lefse", 
                     group = "group", 
                     alpha = 0.05, 
                     lefse_subgroup = NULL,
                     p_adjust_method = "none")

# Figure 4F
t1$plot_diff_cladogram(use_taxa_num = 200, 
                       use_feature_num = 50, 
                       clade_label_level = 5, 
                       group_order = c("Control", "HFD", "FO"),
                       color = c("skyblue","red", "purple")
) 
# Figure 4G
t1$plot_diff_bar(use_number = 1:40, 
                 width = 0.8, 
                 group_order = c("Control", "FO", "HFD"),
                 color_values = c("skyblue", "purple","red"))
