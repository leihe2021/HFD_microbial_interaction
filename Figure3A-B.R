#----------------------------------------------------------------------------------#
# Figure 3A and B 
# Roseburia-colocalized genes are enriched in fatty acid metabolic pathway
#----------------------------------------------------------------------------------#
rm(list = ls())
#------------#
# figure 3A
#------------#
library(tidyverse)
library(clusterProfiler)
library(ggrepel)
input <- read_rds("fig3A_input.rds")

norm <- function(x){
  (x-min(x))/max(x)
}

tmp <- input %>% 
  dplyr::filter(PP4 > 0.3) %>% 
  arrange(chr, start) %>% 
  mutate(
    PP4 = ifelse(exp == "rumi", PP4*-1, PP4),
  ) %>% 
  group_by(chr) %>% 
  mutate(
    start = chr + norm(start),
    start = ifelse(start > 22, 22, start)
  ) %>% 
  ungroup()

ggplot(tmp, aes(x = start, y = PP4, color = exp, fill = exp, shape = exp)) +
  geom_jitter(size = 3) +
  scale_shape_manual(values = c(24,25)) +
  scale_x_continuous(breaks = seq(1,22)) +
  scale_y_continuous(breaks = c(-1.2,-0.3,0.3,1.2)) +
  geom_hline(yintercept = 0, color = "grey", size = 1) +
  ylim(-0.8,1) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Chromosome") +
  geom_text_repel(
    aes(label = symbol),
    size = 3.5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.35, "lines"),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 20))

#----------------------------------------#
# figure 3B colocalized Genes annotation
#----------------------------------------#
kegg.list <- read_rds("colocGenes_KEGG.rds")
dat <- kegg.list$`BRCA|Roseburia`@result
ind <- dat$Description
res.list <- read_rds(file = "brca_rose_keggResList.rds")
df <- lapply(ind, function(x){
  tmp <- res.list[[x]]
  tmp$kegg_pathway <- rep(x, nrow(tmp))
  tmp <- tmp[,c(3,1,2)]
  return(tmp)
}) %>% 
  do.call(rbind, .)

top <- tmp[tmp$exp == "rose",]$symbol
input <- df[df$SYMBOL %in% top,]
sort_cols <- table(input$kegg_pathway) %>% 
  as.data.frame() %>% 
  mutate(pathway = fct_reorder(Var1, .x = Freq))

sort_rows <- table(input$SYMBOL) %>% 
  as.data.frame() %>% 
  mutate(symbol = fct_reorder(Var1, .x = Freq))

input$kegg_pathway <- factor(input$kegg_pathway, levels = levels(sort_cols$pathway))
input$SYMBOL <- factor(input$SYMBOL, levels = levels(sort_rows$symbol))

ggplot(input, aes(x = SYMBOL, y = kegg_pathway, color = kegg_pathway)) +
  geom_point(shape = 15, size = 8) +
  theme_classic() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, size = 12, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_manual(values = cols) +
  ylab("")
ggsave(filename = "scripts/figs_scripts/results_figs_tables/figure_3B_anno_ColocalGenes.pdf",
       width = 11.5, height = 3.3)
