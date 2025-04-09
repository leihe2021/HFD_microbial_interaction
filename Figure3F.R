#----------------------------------------------------------------------------------#
# Figure 3F Manhattan plot for HADHA-rs137852769
#----------------------------------------------------------------------------------#
library(tidyverse)
library(ggtext)
library(normentR)
library(data.table)
sig = 1e-5
sig1 = 1e-8
# filter a small number of SNPs (10% in this case)
gwas <- fread("rs137852769-C100001112.tsv.gz")

sig_data <- gwas %>% 
  subset(pval < sig)

notsig_data <- gwas %>% 
  subset(pval >= sig) %>% 
  group_by(chrom) %>% 
  sample_frac(0.03)

gwas_data <- bind_rows(sig_data, notsig_data)
gwas_data <- na.omit(gwas_data)
gwas_data <- subset(gwas_data, chrom != "X" & chrom != "Y")

unique(gwas_data$chrom)

gwas_data$chrom <- as.numeric(gwas_data$chrom)

data_cum <- gwas_data %>% 
  group_by(chrom) %>% 
  summarise(max_bp = max(pos)) %>% 
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chrom, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chrom") %>% 
  mutate(bp_cum = pos + bp_add)

axis_set <- gwas_data %>% 
  group_by(chrom) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(pval == min(pval)) %>% 
  mutate(ylim = abs(floor(log10(pval))) + 2) %>% 
  pull(ylim)

gwas_data$size <- ifelse(-log10(gwas_data$pval) > 8, "small", "large") 

manhplot <- ggplot(gwas_data, aes(
  x = bp_cum, y = -log10(pval),
  color = as_factor(chrom)
)) +
  geom_hline(
    yintercept = -log10(sig1), color = "grey40",
    linetype = "dashed"
  ) +
  geom_point(aes(size = size), alpha = 0.5) +
  scale_x_continuous(
    label = axis_set$chrom,
    breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(
    c("#7878BA", "#000042"),
    unique(length(axis_set$chrom))
  )) +
  scale_size_manual(values = c(1, 2)) +
  labs(
    x = NULL,
    y = "-log<sub>10</sub>(p)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(size = 8, vjust = 0.5)
  )

