#---------------fig4B  ---------------------#
library(tidyverse)
dat <- readxl::read_xlsx("HFD_model_phenotypedata.xlsx", sheet = 1)
dat[is.na(dat)] <- 0
input <- dat %>% 
  group_by(group) %>% 
  summarise(
    across(.cols = starts_with("day"), list(mean = mean, sd = sd))
  ) %>% 
  gather(key = "key", value = "value", -group) %>% 
  separate(col = key, into = c("remove", "day", "stat_type")) %>% 
  select(-remove) %>% 
  spread(key = stat_type, value = "value") %>% 
  mutate(
    day = as.numeric(day),
    group = factor(group, levels = c("Control", "HFD", "FO"))
  ) %>% 
  arrange(group, day)

input %>% 
  ggplot(aes(x= day, y = mean, color = group, group = group)) +
  geom_line(position = position_dodge(0.2)) +
  geom_point(size = 2.8, shape = 21, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 1.5, 
                position = position_dodge(0.2)) +
  ylab("Tumor size") +
  xlab("Day") +
  #xlim(0, 30) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30)) +
  theme_classic() +
  coord_fixed(ratio = 1/120) +
  scale_color_manual(values = c("#0072B2", "#FF0000", "#B848FF"))

# fig 4B metastasis
library(ggprism)
dat <- readxl::read_xlsx("HFD_model_phenotypedata.xlsx",sheet = 3)
input <- dat %>% 
  select(-no) %>% 
  rename(
    "Liver" = "liver_dis",
    "Lung" = "lung_dis"
  ) %>% 
  mutate(
    group = factor(group, levels = c("Control", "HFD", "FO")) 
  ) %>% 
  gather(key = "metastasis", value = "count", -group)

cols <- c("#0072B2","#FF0000", "#B848FF")

ggplot(input, aes(x = metastasis, y = count, color = group)) +
  geom_boxplot() + 
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(panel.grid = element_blank()) +
  coord_fixed(ratio = 1/8)
