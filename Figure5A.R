#----------------------------
# # Figure 5 A
# metabolism data
#----------------------------
library(tidyverse)
library(scatterplot3d)
dat <- read.csv("PCA_with_All_Samples_Score.csv")
input <- dat %>% 
  filter(!grepl("QC", SampleID)) %>% 
  mutate(
    group = case_when(
      grepl("L\\d+", SampleID) ~ "Con",
      grepl("F\\d+", SampleID) ~ "FOD",
      grepl("H\\d+", SampleID) ~ "HFD"
    )
      ) %>% 
  mutate(group = factor(group, levels = c("Con", "HFD", "FOD"))) %>% 
  select(group, everything()) %>% 
  arrange(group)

# Figure 5A 3D PCAplot for metabolites  
color = c("skyblue","red", "purple")
cols <- color[as.numeric(input$group)]

pdf(file = "metabolism_pca.pdf", width = 6.2, height = 6)
scatterplot3d(input[,3:5], pch = 16, color = cols, size= 10, cex.symbols = 1.5)
legend(
  "bottom", legend = c("Con","HFD","FOD"),
  col =  cols, 
  pch = c(16, 16, 16), 
  inset = -0.25, xpd = TRUE, horiz = TRUE)
dev.off()

