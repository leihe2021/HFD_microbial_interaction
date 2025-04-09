#----------------------------
#' metabolism data
#' Figure 5B
#' BRCA mouse model data
#----------------------------
library(tidyverse)
ind <- c("CON_vs_HFD", "FO_vs_HFD", "FO_vs_CON")
res <- lapply(ind, function(x){
  dir = paste("Treatment/", AllMet_Tests_Joined.csv", sep = x)
  tmp <- read_csv(file = dir)
  dat <- tmp %>% 
    select(Metabolite, log2FC, Uni_P, Uni_FDR, OPLSDA_VIP) %>% 
    mutate(
      sig = case_when(
        Uni_P < 0.05 & abs(log2FC) > 0.53 & OPLSDA_VIP > 1 ~ 1,
        TRUE ~ 0
      ),
      grp = x
    ) %>% 
    select(grp, everything())
})
names(res) <- ind

# plot combined volcano
library(tidyverse)
library(ggrepel)
dat.list <- read_rds("plot_InputData/fig6_volcanoInput.rds")

dat <- lapply(1:3, function(x){
  df <- dat.list[[x]] %>% 
    dplyr::rename(pval = Uni_P) %>% 
    mutate(
      significance = case_when(
        log2FC > 0.53 & log2FC < 1 & pval < 0.05 & OPLSDA_VIP > 1 ~ "up",
        log2FC < -0.53 & log2FC > -1 & pval < 0.05 & OPLSDA_VIP > 1 ~ "down",
        log2FC > 1 & pval < 0.05 & OPLSDA_VIP > 1 ~ "up_ad",
        log2FC < -1 & pval < 0.05 & OPLSDA_VIP > 1 ~ "down_ad",
        TRUE ~ "insig"
      ),
      x.local = x,
      grp = factor(grp, levels = c("CON_vs_HFD", "FO_vs_HFD", "FO_vs_CON")),
      pval_new = log2(pval)/10 + 3*x - 1,
      significance = factor(significance, levels = c("up_ad", "up","insig", "down", "down_ad"))
    ) %>% 
    select(Metabolite, log2FC, pval, pval_new, grp, significance, x.local)
  return(df)
}) %>% 
  do.call(rbind, .)

x_local_down <- dat %>% 
  group_split(grp) %>% 
  lapply(function(x){range(x$pval)[1]})

x_local_up <- dat %>% 
  group_split(grp) %>% 
  lapply(function(x){range(x$pval)[2]})


bg <- data.frame(grp = levels(dat$grp),x_down = c(0,3,6), x_up = c(2,5,8),
                 y_down = -6, y_up = 6)

dat.label <- dat %>% 
  filter(significance %in% c("up_ad", "down_ad")) %>% 
  mutate(
    sig = case_when(
      log2FC > 2.3 ~ "up",
      log2FC < -2.3 ~ "down",
      TRUE ~ "no"
    )
  ) %>% 
  filter(sig != "no")


p1 <- ggplot() +
  geom_jitter(dat, mapping = aes(x = pval_new, y = log2FC, color = significance), 
              width = 0.3) +
  scale_color_manual(values = c("red","#FAC0AE","grey80","skyblue", "#6175DB")) +
  annotate("rect", xmin = 0, xmax = 2.5, ymin = -6, ymax = 6, alpha = .1, fill = "grey50") +
  annotate("rect", xmin = 3, xmax = 5.5, ymin = -6, ymax = 6, alpha = .1, fill = "grey50") +
  annotate("rect", xmin = 6, xmax = 8.5, ymin = -6, ymax = 6, alpha = .1, fill = "grey50") +
  annotate("rect", xmin = 0, xmax = 2.5, ymin = -6.5, ymax = -5.5, alpha = .8, fill = "#B95D43") +
  annotate("rect", xmin = 3, xmax = 5.5, ymin = -6.5, ymax = -5.5, alpha = .8, fill = "#F43981") +
  annotate("rect", xmin = 6, xmax = 8.5, ymin = -6.5, ymax = -5.5, alpha = .8, fill = "#7938C3") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) 
p1

p1 + geom_text_repel(
  data = dat.label, aes(x = pval_new, y =  log2FC, label = Metabolite),
  size = 3,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.35, "lines"),
  max.overlaps = getOption("ggrepel.max.overlaps", default = 20)
)

