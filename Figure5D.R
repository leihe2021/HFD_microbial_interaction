#------------------------------------------------------
# # Figure 5D microbe-metabolites correlation analysis
# metabolism data
#------------------------------------------------------
rm(list = ls())
library(ggprism)
library(ggpubr)
library(tidyverse)
library(cowplot)
load("metabolism_annlysisData.Rda")
met <- m_met %>% 
  column_to_rownames(var = "Metabolite") %>% 
  t() %>% 
  as.data.frame()

input <- met %>% 
  rownames_to_column(var = "group") %>% 
  as.data.frame() %>% 
  mutate(
    group = case_when(
      grepl("L", group) ~ "Con",
      grepl("H", group) ~ "HFD",
      grepl("F", group) ~ "FOD"
    )
  ) %>% 
  gather(key = "metabolites", value = "acc", -group) %>% 
  select(group, everything()) %>% 
  mutate(
    group = factor(group, levels = c("Con", "HFD", "FOD"))
  ) %>% 
  arrange(metabolites, group)

# figure 5D 
barplot_metabolites <- function(dat, met){
  dat %>% 
    filter(metabolites == met) %>% 
    ggplot(aes(x = group, y = acc, color = group)) +
    geom_boxplot() +
    theme_classic() +
    scale_color_manual(values = color) +
    theme(legend.position = "none") +
    ylab(met) +
    xlab("") 
}
color = c("skyblue","red", "purple")  
ind <- c("Leucine", "Isoleucine", "Valine", "Ketoleucine")
plist <- lapply(ind, function(x)barplot_metabolites(dat = input, met = x))
names(plist) <- ind
cowplot::plot_grid(plotlist = plist, ncol = 2)
