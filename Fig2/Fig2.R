library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(here)

measurements <- read_csv(here("Fig2", "measurements.csv")) %>%
  janitor::clean_names()

measurements$tissue <- factor(
  measurements$tissue,
  levels = c("Liver", "Intestine"),
  # labels = c("Liver Miracidia", "Intestine Miracidia")
)

#---------------Area---------------

df_area <- ggplot(measurements, aes(x = tissue, y = area_1_blue_um2, fill = tissue)) +
  geom_boxplot(alpha = 1.0, size = 1.0, outliers = FALSE) +
  geom_quasirandom(color = "black", size = 2, alpha = 0.6) +
  stat_summary(geom = 'point', fun = mean, size = 3, color = 'grey', shape = 'triangle') +
  theme_minimal() +
  scale_fill_manual(values = c("Liver" = "steelblue", "Intestine" = "indianred")) +
  stat_compare_means(method = 't.test', size = 3, label.x.npc = .5, hjust = 0.5) +
  labs(x = "Tissue source", y = "Miracidia area (µm²)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
df_area

#---------------Length---------------

df_length <- ggplot(measurements, aes(x = tissue, y = length_1_blue_um, fill = tissue)) +
  geom_boxplot(alpha = 1.0, size = 1.0, outliers = FALSE) +
  geom_quasirandom(color = "black", size = 2, alpha = 0.6) +
  stat_summary(geom = 'point', fun = mean, size = 3, color = 'grey', shape = 'triangle') +
  theme_minimal() +
  scale_fill_manual(values = c("Liver" = "steelblue", "Intestine" = "indianred")) +
  stat_compare_means(method = 't.test', size = 3, label.x.npc = .5, hjust = 0.5) +
  labs(x = "Tissue source", y = "Miracidia length (µm)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

combined_plot <- plot_grid(df_area, df_length, labels = c("A", "B"), ncol = 2)
combined_plot

save_plot(here("Fig2", "Fig2.pdf"), combined_plot, base_height = 4, base_width = 6.5)
save_plot(here("Fig2", "Fig2.png"), combined_plot, base_height = 4, base_width = 6.5, bg = 'white')

#---------------Summary---------------

summary_table <- measurements %>%
  group_by(tissue) %>%
  summarise(
    area_mean = mean(area_1_blue_um2, na.rm = TRUE),
    length_mean = mean(length_1_blue_um, na.rm = TRUE),
  )
