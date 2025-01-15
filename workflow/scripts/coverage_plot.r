library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(caret)
library(viridis)

covg_plot <- read.csv('/nas/longleaf/home/adaigle/TEforest/workflow/scripts/covg_exp.csv')

covg_plot_long <- covg_plot %>%
  pivot_longer(cols = c(F1, Precision, Recall), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, "F1" = "F1 Score"))

colors <- c(
  "TEforest" = "#F8766D", 
  "temp2" = "#CD9600", 
  "retroseq" = "#7CAE00", 
  "temp" = "#00BE67", 
  "popoolationte" = "#00BFC4", 
  "teflon" = "#00A9FF", 
  "popoolationte2" = "#C77CFF", 
  "tepid" = "#FF61CC"
)
# Reorder the levels of the caller factor to match the order in colors
covg_plot_long$caller <- factor(covg_plot_long$caller, levels = names(colors))

# Create the plot
cov <- ggplot(covg_plot_long, aes(x = Coverage, y = value, color = caller)) +
  geom_line(size = 1.5, alpha = 0.5) +
  geom_point(size=3.5) +
  facet_wrap(~ metric) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50), labels = c(5, 10, 20, 30, 40, 50)) +
  labs(title = "",
       x = "Coverage",
       y = "Value") +
  theme_minimal() + 
  scale_color_manual(values = colors) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 18), 
        strip.text = element_text(size = 18))
ggsave("/nas/longleaf/home/adaigle/TEforest/plots/cov_150.png", plot = cov, width = 12, height = 6, dpi = 300)

