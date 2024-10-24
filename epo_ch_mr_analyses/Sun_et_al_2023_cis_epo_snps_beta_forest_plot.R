### Set working directory and load required packages:
setwd()
library(tidyverse)
library(vroom)

### Load data:
epo_snps <- vroom("epo_snps.csv")

### Plot for epo_hgb_chr7:
epo_snps$variant <- factor(epo_snps$variant, levels = c("rare SNP: rs11976235-T", "common SNP: rs1617640-A"))
plot_epo_snps <- ggplot(data = epo_snps, aes(x = SNP, y = B, ymin = lower_CI, ymax = upper_CI)) +
  geom_pointrange(aes(col = variant), shape = 15, size = 1.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(NULL) +
  ylab("Beta (95% CI)") +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI, col = variant), width = 0.3, cex = 2.00) + 
  facet_wrap(~phenotype, strip.position="left", nrow=9, scales = "free_y") +
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    strip.text.y = element_text(hjust=0.5, vjust = 0.5, angle = 90),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(color = "black", fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 20),
    legend.background = element_rect(color = "black", size = 1.0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
  ) + 
  coord_flip() +
  guides(color = guide_legend(reverse = TRUE))
plot_epo_snps

### Save plot:
ggsave(plot_epo_snps, file = "Sun_et_al_2023_cis_epo_snps_forest_plot.pdf")
