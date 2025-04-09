#-----------------------------------------------------------------------#
# Figure 2 outside barplot 24 representative UVMR results               
#-----------------------------------------------------------------------#
library(tidyverse)
library(cowplot)

uvmr_list <- read_rds("genus_MR_result.rds")

uvmr <- lapply(names(uvmr_list), function(x){
  tmp <- uvmr_list[[x]] %>% 
    filter(pval < 0.05) %>% 
    arrange(exposure) %>% 
    mutate(type = x,
           trait = paste(cancer_type, exposure, sep = "|")) %>% 
    select(type, exposure, trait, pval, everything())
}) %>% 
  do.call(rbind, .)

ind <- as.data.frame(table(uvmr$trait)) %>% 
  rename(trait = Var1) %>%
  arrange(desc(Freq)) %>% 
  filter(Freq > 1) %>% 
  mutate(trait = as.character(trait)) %>% 
  pull(trait)


df <- do.call(rbind, uvmr_list) %>% 
  mutate(trait = paste(cancer_type, exposure, sep = "|")) %>% 
  select(trait, or, pval, method) %>% 
  filter(trait %in% ind) %>% 
  mutate(lgOR  = log2(or)) %>% 
  mutate(sig = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 & pval >= 0.001 ~ "**",
    pval < 0.05 & pval >= 0.01 ~ "*",
    T ~ ""
  ))

rownames(df) <- NULL
ind <- sort(ind)

p.list <- lapply(ind, function(x){
  df %>% 
    filter(trait == x) %>% 
    ggplot(aes(x = method, y = or)) +
    geom_col(aes(fill = method), width = 0.8) +
    scale_fill_manual(values = c("#0078E9", "#FFA500", "purple", "#C40000", "#700045")) +
    geom_hline(yintercept = 1, color = "white", linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_text(aes(y = or + 0.02, label = sig), size = 8) +
    ylab(paste(gsub("(.*)\\.+(.*?)", "\\2", x), "(OR)", sep = " ")) +
    ggtitle(label = gsub("(.*?)\\|(.*)", "\\1", x)) +
    xlab("")
})

p.list[[1]]

plot_grid(plotlist = p.list, ncol = 6)
ggsave(filename = "figure_2_outside_barplot.pdf" 
       ,width = 15, height = 10)

#-------------------------------------------------------------#
# figure 2 inner circos plot 
#-------------------------------------------------------------#
# circos input data preparing
# circos plot ")
library(tidyverse)
library(TwoSampleMR)
# load the result data of MVMR
mvmr <- read_rds("clean_merge_MVMResult.rds")

mv_genus <- sapply(strsplit(as.character(mvmr$trait), split = "\\|"), "[", 2)
mv_genus <- mv_genus[!duplicated(mv_genus)]

# selected high frequent genus to plot circos
uvmr_total_filtered <- read_rds("uvmr_results.rds")

uvmr <- readxl::read_xlsx("genus_mr_p005_results.xlsx") %>% 
  select(cancer_type, outcome, exposure, method, or, or_lci95, or_uci95, pval) %>% 
  as.data.frame() %>% 
  group_by(cancer_type, exposure) %>% 
  arrange(cancer_type, exposure, pval) %>% 
  slice(1)
dt <- as.data.frame(table(uvmr$exposure)) %>% 
  filter(Freq > 5) %>% 
  rename(genus = Var1, freq = Freq)
dt$genus <- gsub("genus\\.+", "", dt$genus)
genus <- c(dt$genus, mv_genus)
genus <- genus[!duplicated(genus)]
gen <- c(mv_genus, setdiff(genus, mv_genus))
input <- uvmr_total_filtered %>% 
  mutate(
    genus = gsub("genus\\.+", "", exposure)
  ) %>% 
  filter(genus %in% gen) %>% 
  mutate(genus = factor(genus, levels = gen)) %>% 
  select(outcome, cancer_type, exposure, genus, everything())

input$sig <- ifelse(input$pval < 0.05, T, F)  
sum(input$sig)  

# set input files
all_names <- names(table(input$cancer_type))
len_res <- as.vector(table(input$cancer_type))
all_species <- gen

# file 1 construct mega_data
df1 <- data.frame(
  v1 = "chr",
  v2 = "-",
  vi = paste0("c", seq_along(all_names)),
  v3 = all_names,
  v4 = 1,
  v5 = len_res,
  v6 = paste0("chrom", seq_along(all_names))
) %>% 
  arrange(desc(row_number()))

df2 <- data.frame(
  v1 = "chr",
  v2 = "-",
  vi = paste0("s", seq_along(all_species)),
  v3 = all_species, 
  v4 = 1, 
  v5 = 2, 
  v6 = paste0("species", seq_along(all_species))
)
df <- rbind(df1,df2)

write.table(df, file = "circos_mega_data.txt", quote = F, 
            col.names = F, row.names = F, sep = " ")

# file 2 set color
color.ls <- c("#D2001A", "blue", "#5800FF", "#FFB200",
              "#FA2FB5", "skyblue","#7DCE13", "#8758FF",
              "#66BFBF", "#FE2A2A", "#FEA704", "#2747B3",
              "#289A48", "#7D4D3C", "#6C286C", "grey")

names(color.ls) <- c(all_names, "not_sig")
scales::show_col(color.ls)
color_df <- col2rgb(color.ls)
color_table <- c()
for (i in 1:(ncol(color_df)-1)) {
  tmp <- paste0("chrom",
                i,
                "=rgb(",
                color_df[1,i],
                ",",
                color_df[2,i],
                ",",
                color_df[3,i],
                ")"
  )
  color_table <- c(color_table, tmp)
}
color_table

writeLines(color_table, "circos_color_table.txt")

# file 3 set links
name_res <- input %>% 
  as.data.frame() %>% 
  select(cancer_type, genus) %>% 
  group_split(cancer_type)

ind <- lapply(name_res, function(x)unique(x$cancer_type)) %>% 
  unlist

name_res <- lapply(name_res, function(x)x$genus)
names(name_res) <- ind

link_df <- data.frame()

for(i in 1:length(all_names)){
  this_name <- paste0("c",i)
  this_cancer <- df %>% 
    filter(vi == this_name) %>% 
    pull(v3)
  
  this_cancer_species <- name_res[[this_cancer]]
  this_cancer_species_idx <- df2 %>% 
    filter(v3 %in% this_cancer_species) %>% 
    pull(vi)
  
  this_link <- data.frame(
    v1 = this_name,
    v2 = seq_along(this_cancer_species_idx) - 1,
    v3 = seq_along(this_cancer_species_idx),
    v4 = rev(this_cancer_species_idx),
    v5 = 1,
    v6 = 2
  )
  link_df = rbind(link_df, this_link)
}
write.table(link_df, file = "circos_mega_link.txt", 
            quote = F, col.names = F,
            row.names = F, sep = " ")

# file 4 set species label
#label_df <- data.frame()
label_df <- df2 %>% 
  select(vi, v4, v5, v3)

#head(label_df)
write.table(label_df, file = "circos_mega_label.txt",
            quote = F, col.names = F, row.names = F,
            sep = " ")

#----------------------------#
# figure 2 inner perl circos
#----------------------------#
# copy above files in dir fold
# cd circos/
# perl "circos-0.69/bin/circos" -config  "../fig2_circos/circos.config"






