# Figure 4H
# BRCA mouse model data
# Co-occurrence analysis
library(tidyverse)
library(igraph)
library(magrittr)
library(plyr)

load("genus_groups_correlation.Rda")
load("otu_grp_taxa.Rda")
genus <- read_rds("genus_relative_accumulation.rds")
tmp <- res.list$con %>% 
  separate(col = genus_pair, into = c("from", "to"), sep = "_") %>% 
  filter(from == "Roseburia" | to == "Roseburia") %>% 
  mutate(sig = abs(correlation) > 0.6) %>% 
  filter(sig)

ind <- unique(c(tmp$from, tmp$to))
genus <- genus[,ind]

get_net <- function(genus_accu, grp = "LFD", res.list, taxa, sel_genus_index){
  links <- res.list[[grp]] %>% 
    separate(col = genus_pair, into = c("from", "to"), sep = "_") %>% 
    filter(from %in% sel_genus_index, to %in% sel_genus_index) %>% 
    mutate(weight = -log2(pval)) %>% 
    mutate(cls = case_when(
      correlation > 0.9 ~ "#CB0000",
      correlation >= 0.6 & correlation <= 0.9 ~ "#FE8B6A",
      correlation < 0.6 & correlation > -0.6 ~ NA,
      correlation <= -0.6 & correlation >= -0.9 ~ "#C1D8FF",
      correlation < -0.9 ~ "#A9CBFF"
    ))
  nodes <- taxa %>% 
    mutate(g = gsub("g__", "", genus)) %>% 
    select(phylum, g) %>% 
    filter(g %in% sel_genus_index) %>% 
    distinct(phylum, g) %>% 
    dplyr::rename(media = g, type_label = phylum) %>% 
    select(media, type_label) %>% 
    mutate(type_label = gsub("p__", "", type_label)) %>% 
    left_join(tmp, by = "media") %>% 
    mutate(
      id = paste("S", 1:length(sel_genus_index), sep = "")
    ) %>% 
    select(id, everything()) 
  links$from <- mapvalues(links$from, nodes$media, nodes$id)
  links$to <- mapvalues(links$to, nodes$media, nodes$id)  
   net <- graph_from_data_frame(d = links, vertices = nodes, directed = T) 
  res = list(net = net, links = links, nodes = nodes)
  return(res)
}
res <- get_net(genus_accu = genus, grp = "Con", res.list = res.list, taxa = taxa, sel_genus_index = ind)
res <- get_net(genus_accu = genus, grp = "HFD", res.list = res.list, taxa = taxa, sel_genus_index = ind)
res <- get_net(genus_accu = genus, grp = "FO", res.list = res.list, taxa = taxa, sel_genus_index = ind)

net <- res$net
links <- res$links
nodes <- res$nodes
V(net)$color <- revalue(nodes$type_label, c("Firmicutes" = "#0078E9",
                                            "Bacteroidota" = "#C40000",
                                            "Proteobacteria" = "#FFA500"))
V(net)$size <- V(net)$size *3
V(net)$label <- V(net)$media
E(net)$arrow.mode <- 0
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1 + E(net)$weight/6 

pdf(file = "cooccurence_3.pdf", width = 8, height = 8)
plot(net,
     layout = layout_in_circle, 
     edge.curved = .1, 
     vertex.label.color = V(net)$color, 
     vertex.label.dist = -2, 
     edge.color = links$cls)
dev.off()
