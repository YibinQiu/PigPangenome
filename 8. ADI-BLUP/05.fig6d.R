rm(list = ls())
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(RColorBrewer)

data <- read_excel("ADI-BLUP.xlsx",sheet = 3) 
data1 <- melt(data, id.vars = c("Trait","method"), 
              variable.name = "Model", value.name = "Accuracy")

data1$interaction <- paste(data1$Model, data1$method, sep = ":")

data1$interaction <- factor(data1$interaction, 
                            levels = c("GBLUP:50K", "ADI-BLUP:50K", 
                                       "GBLUP:PanSNP", "ADI-BLUP:PanSNP", 
                                       "GBLUP:PanSV", "ADI-BLUP:PanSV",
                                       "GBLUP:PanINDEL", "ADI-BLUP:PanINDEL"))
data1 <- data1 %>%
  group_by(Trait, method) %>%
  mutate(group_line = paste(Trait, method))

ggplot(data1, aes(x = interaction, y = Accuracy)) +
  geom_boxplot(aes(color = Model), outlier.size = 0, linewidth = 0.8, width = 0.5, fill = NA) +  
  geom_point(aes(color = Model), size = 1.5) +  
  geom_line(aes(group = group_line), color = "gray70", linewidth = 0.2, na.rm = TRUE) +  
  labs(x = NULL, y = "Prediction Accuracy") +
  scale_y_continuous(breaks = seq(0.1, 0.6, by = 0.1)) +
  scale_color_manual(values = c("#e6846D", "#8DCDD5")) +  
  theme_minimal() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(color = "black", angle = -45, hjust = 0), 
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())
ggsave("Acc.pdf", width = 8,height = 6)
